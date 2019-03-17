Shader "Unlit/Analytical Normal"
{
    Properties
    {
        _MainTex ("Texture", 2D) = "white" {}
    }
    SubShader
    {
        Tags { "RenderType"="Opaque" }
        LOD 100

        Pass
        {
            CGPROGRAM
            #pragma vertex vert
            #pragma fragment frag
            // make fog work
            #pragma multi_compile_fog

            #include "UnityCG.cginc"

            struct appdata
            {
                float4 vertex : POSITION;
                float2 uv : TEXCOORD0;
            };

            struct v2f
            {
                float2 uv : TEXCOORD0;
                UNITY_FOG_COORDS(1)
                float4 vertex : SV_POSITION;
            };

            sampler2D _MainTex;
            float4 _MainTex_ST;

            v2f vert (appdata v)
            {
                v2f o;
                o.vertex = UnityObjectToClipPos(v.vertex);
                o.uv = TRANSFORM_TEX(v.uv, _MainTex);
                UNITY_TRANSFER_FOG(o,o.vertex);
                return o;
            }

			float hash(float n) { return frac(sin(n)*753.5453123); }


			//---------------------------------------------------------------
			// value noise, and its analytical derivatives
			//---------------------------------------------------------------

			float4 noised(in float3 x)
			{
				float3 p = floor(x);
				float3 w = frac(x);
				float3 u = w * w*(3.0 - 2.0*w);
				float3 du = 6.0*w*(1.0 - w);

				float n = p.x + p.y*157.0 + 113.0*p.z;

				float a = hash(n + 0.0);
				float b = hash(n + 1.0);
				float c = hash(n + 157.0);
				float d = hash(n + 158.0);
				float e = hash(n + 113.0);
				float f = hash(n + 114.0);
				float g = hash(n + 270.0);
				float h = hash(n + 271.0);

				float k0 = a;
				float k1 = b - a;
				float k2 = c - a;
				float k3 = e - a;
				float k4 = a - b - c + d;
				float k5 = a - c - e + g;
				float k6 = a - b - e + f;
				float k7 = -a + b + c - d + e - f - g + h;

				return float4(k0 + k1 * u.x + k2 * u.y + k3 * u.z + k4 * u.x*u.y + k5 * u.y*u.z + k6 * u.z*u.x + k7 * u.x*u.y*u.z,
					du * (float3(k1, k2, k3) + u.yzx*float3(k4, k5, k6) + u.zxy*float3(k6, k4, k5) + k7 * u.yzx*u.zxy));
			}

			float4 sdBox(float3 p, float3 b) // distance and normal
			{
				float3 d = abs(p) - b;
				float x = min(max(d.x, max(d.y, d.z)), 0.0) + length(max(d, 0.0));
				float3  n = step(d.yzx, d.xyz)*step(d.zxy, d.xyz)*sign(p);
				return float4(x, n);
			}

			float4 fbmd(in float3 x)
			{
				const float scale = 1.5;

				float a = 0.0;
				float b = 0.5;
				float f = 1.0;
				float3  d = float3(0.0, 0.0, 0.0);
				for (int i = 0; i < 8; i++)
				{
					float4 n = noised(f*x*scale);
					a += b * n.x;           // accumulate values		
					d += b * n.yzw*f*scale; // accumulate derivatives
					b *= 0.4;             // amplitude decrease
					f *= 1.8*cos(_Time.y * 0.25);             // frequency increase
				}

				return float4(a, d);
			}

			float4 map(in float3 p)
			{
				float4 d1 = fbmd(p);
				//d1.x -= 0.37;
				d1.x -= 0.37;
				d1.x *= 0.7;
				d1.yzw = normalize(d1.yzw);

				// clip to box
				float4 d2 = sdBox(p, float3(1.5, 1.5, 1.5));
				return (d1.x > d2.x) ? d1 : d2;
			}

			float2 iBox(in float3 ro, in float3 rd, in float3 rad)
			{
				float3 m = 1.0 / rd;
				float3 n = m * ro;
				float3 k = abs(m)*rad;
				float3 t1 = -n - k;
				float3 t2 = -n + k;
				float tN = max(max(t1.x, t1.y), t1.z);
				float tF = min(min(t2.x, t2.y), t2.z);
				if (tN > tF || tF < 0.0) return float2(-1.0, -1.0);
				return float2(tN, tF);
			}

			// raymarch
			float4 interesect(in float3 ro, in float3 rd)
			{
				float4 res = float4(-1.0, -1.0, -1.0, -1.0);

				// bounding volume    
				float2 dis = iBox(ro, rd, float3(1.5, 1.5, 1.5));
				if (dis.y < 0.0) return res;

				// raymarch
				float tmax = dis.y;
				float t = dis.x;
				for (int i = 0; i < 128; i++)
				{
					float3 pos = ro + t * rd;
					float4 hnor = map(pos);
					res = float4(t, hnor.yzw);

					if (hnor.x < 0.001) break;
					t += hnor.x;
					if (t > tmax) break;
				}

				if (t > tmax) res = float4(-1.0, -1.0, -1.0, -1.0);
				return res;
			}

			// fibonazzi points in s aphsre, more info:
			// http://lgdv.cs.fau.de/uploads/publications/spherical_fibonacci_mapping_opt.pdf
			float3 forwardSF(float i, float n)
			{
				static const float PI = 3.141592653589793238;
				static const float PHI = 1.618033988749894848;
				float phi = 2.0*PI*frac(i / PHI);
				float zi = 1.0 - (2.0*i + 1.0) / n;
				float sinTheta = sqrt(1.0 - zi * zi);
				return float3(cos(phi)*sinTheta, sin(phi)*sinTheta, zi);
			}

			float calcAO(in float3 pos, in float3 nor)
			{
				float ao = 0.0;
				for (int i = 0; i < 32; i++)
				{
					float3 ap = forwardSF(float(i), 32.0);
					float h = hash(float(i));
					ap *= sign(dot(ap, nor)) * h*0.25;
					ao += clamp(map(pos + nor * 0.001 + ap).x*3.0, 0.0, 1.0);
				}
				ao /= 32.0;

				return clamp(ao*5.0, 0.0, 1.0);
			}

            fixed4 frag (v2f i) : SV_Target
            {
				/*
                // sample the texture
                fixed4 col = tex2D(_MainTex, i.uv);
                // apply fog
                UNITY_APPLY_FOG(i.fogCoord, col);
                return col;
				*/

				float2 p = (-_ScreenParams.xy + 2.0*(_ScreenParams.xy * i.uv)) / _ScreenParams.y;

				// camera anim
				float an = sin(0.98803162);// 0.3*_Time.y;
				//float an = 0.3*_Time.y;
				float3 ro = 4.5*float3(cos(an), 0.5, sin(an));
				float3 ta = float3(0.0, 0.0, 0.0);

				// camera matrix	
				float3  cw = normalize(ta - ro);
				float3  cu = normalize(cross(cw, float3(0.0, 4.5, 0.0)));
				float3  cv = normalize(cross(cu, cw));
				//float3  rd = normalize(p.x*cu + p.y*cv + 1.7*cw);
				float3  rd = normalize(p.x*cu + p.y*cv + 1.7*cw);

				// render
				float3 col = float3(0.0, 0.0, 0.0);
				float4 tnor = interesect(ro, rd);
				float t = tnor.x;

				if (t > 0.0)
				{
					float3 pos = ro + t * rd;
					#ifndef SHOW_NUMERICAL_NORMALS
					float3 nor = tnor.yzw; // no need to call calcNormal( pos );
					#else
					float3 nor = calcNormal(pos);
					#endif
					float occ = calcAO(pos, nor);
					float fre = clamp(1.0 + dot(rd, nor), 0.0, 1.0);
					float fro = clamp(dot(nor, -rd), 0.0, 1.0);
					//col = lerp(float3(0.05, 0.2, 0.3), float3(1.0, 0.95, 0.85), 0.5 + 0.5*nor.y);
					col = lerp(float3(0.05, 0.2, 0.3), float3(0.95, 0.95, 0.95), 0.5 + 0.5*nor.y);
					//col = 0.5+0.5*nor;
					col += 10.0*pow(fro, 12.0)*(0.04 + 0.96*pow(fre, 5.0));
					col *= pow(float3(occ, occ, occ), float3(1.0, 1.1, 1.1));
				}

				col = sqrt(col);
				return float4(col, 1.0);
				//fragColor = vec4(col, 1.0);
            }
            ENDCG
        }
    }
}
