const vec3 CAMERA_POS = vec3(0, 2.5, -10);

const vec3 LIGHT_POS_1 = vec3(4, 0, 0);
const float LIGHT_RAD_1 = 0.2;
const vec3 LIGHT_COL_1 = vec3(1,1,1);
const float LIGHT_POW_1 = 2.0;

const vec3 LIGHT_POS_2 = vec3(-3, 0.8, 2.5);
const float LIGHT_RAD_2 = 0.2;
const vec3 LIGHT_COL_2 = vec3(1,1,1);
const float LIGHT_POW_2 = 10.0;

const vec3 LIGHT_POS_3 = vec3(0, 1.5, 0);
const float LIGHT_RAD_3 = 0.0;
const vec3 LIGHT_COL_3 = vec3(1,0,0);
const float LIGHT_POW_3 = 9.0;

const float INF = 1e10;
const float EPS = 1e-6;
const float AIR_N = 1.0003;
const float GLASS_N = 1.51;
const float DIAMOND_N = 2.5;

const float K = 0.1;
const float HASH_VAL =  0.1031;
vec3 randDir;

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

vec3 fade(vec3 t)
{ 
    return t*t*t*(t*(6.*t-15.)+10.); 
}

float hash(vec3 p3)
{
	p3 = fract(p3 *HASH_VAL); 
    p3 += dot(p3, p3.yzx + 19.19);
    return fract((p3.x + p3.y) * p3.z);
}

float grad(float hash, vec3 p) 
{
    int h = int(1e4*hash) & 15;
	float u = h<8 ? p.x : p.y,
 		  v = h<4 ? p.y : h==12||h==14 ? p.x : p.z;
    return ((h&1) == 0 ? u : -u) + ((h&2) == 0 ? v : -v);
}

float powf(float x, int a) 
{
    float q = x;
    int i = a;
    while (1 < i--) 
    {
        q *= x;
    }
    return a == 0 ? 0.0 : q;
}

float noise(vec3 x)            
{
    vec3 p = x;
    vec3 i = floor(p);
	vec4 a = dot(i, vec3(1, 57, 21)) + vec4(0, 57, 21, 78);
	vec3 f = cos((p-i)*acos(-1.0))*(-0.5) + 0.5;
	a = mix(sin(cos(a)*a),sin(cos(1.0 + a)*(1.0 + a)), f.x);
	a.xy = mix(a.xz, a.yw, f.y);
	return mix(a.x, a.y, f.z);
}

float rand(float frame)    
{
    return fract(sin( dot(vec3(frame), vec3(23.4653,67.907,43.89037) )) * 12903.5909);
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

struct MatType 
{
    float EMIS;
    float DIFF;
    float REFL;
    float REFR;
    float N;
    vec3  COL;
};

struct Sphere 
{
    vec3  POS;
    float RAD;
    MatType M;
};

struct LightSrc 
{
    Sphere SPH;
    float  POW;
};

struct Triangle 
{
    vec3 VRT[3];
};

struct Pyramid                   
{
    vec3    VRT[5];
    MatType M;
};

struct Scene 
{
    LightSrc[3] L;
    Pyramid     P;
    Sphere      S;
} SC;

MatType Emissive = MatType(1.0, 0.0, 0.0, 0.0, AIR_N, vec3(0));
MatType Diffusive = MatType(0.0, 1.0, 0.0, 0.0, AIR_N, vec3(0));
MatType Reflective = MatType(0.0, 0.0, 1.0, 0.0, AIR_N, vec3(0));
MatType Refractive = MatType(0.0, 0.0, 0.0, 1.0, GLASS_N, vec3(0));
MatType MyGlass = MatType(0.0, 0.4, 1.0, 1.0, GLASS_N, vec3(0.0,4.0,0.0));

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

float TraceFire(vec3 p)
{
    return length(vec3(0, 1, 0) - p*vec3(1, 0.7, 1)) - 1.0 +
        (noise(p + vec3(0, 5, 0)) + noise(p * 3.0)* 0.5)* 0.25 *(p.y) ;
}

vec4 FireRayMarch(vec3 org, vec3 dir)
{
	float d = 2.1, glow = 1.1, eps = 0.01;
	vec3  p = org;
	bool glowed = false;
    float k = 64.0;
	
	for(int i = 0; i < int(k); i++)
	{
		d = min(150.0 - length(p), abs(TraceFire(p))) + eps;
		p += d * dir;
		if (d > eps)
		{
            glowed = glowed || TraceFire(p) < 0.0;
			glow = glowed ? float(i)/k : 0.0;
		}
	}
	return vec4(p,glow);
}

float TracePlane(vec3 pos, vec3 dir, out vec3 normal)
{
    float t = (-1.2 - pos.y) / dir.y;
    if (t <= 0.0) {
        return INF;
    }
    
    vec3 worldPos = t * dir + pos;
    if (dot(worldPos.xz, worldPos.xz) >= 100.0) {
        return INF;
    }
    normal = vec3(0, 1, 0);
    return t;
}

float TraceTriangle(vec3 pos, vec3 dir, Triangle T, out vec3 N)
{
    vec3 E1 = T.VRT[1] - T.VRT[0];
    vec3 E2 = T.VRT[2] - T.VRT[0];
    vec3 H = cross(dir, E2);
    float a = dot(E1, H);
    float b = 1.0/a;
    vec3 S = pos - T.VRT[0];
    float u = b * dot(S, H);
    vec3 Q = cross(S, E1);
    float v = b * dot(dir, Q);
    float t = b * dot(E2, Q);
    
    if ((a > -EPS && a < EPS) || u < 0.0 || 
        u > 1.0 || v < 0.0 || u + v > 1.0 || t < EPS) 
    {
        return INF;
    }
    
    vec3 V1 = T.VRT[2] - T.VRT[0];
    vec3 V2 = T.VRT[1] - T.VRT[0];
    
    N = normalize(cross(V1,V2));
    
    return t;
}

float TracePyramid(vec3 pos, vec3 dir, Pyramid P, out vec3 N) 
{
    Triangle T1 = Triangle(vec3[3](P.VRT[0], P.VRT[1], P.VRT[2]));
    Triangle T2 = Triangle(vec3[3](P.VRT[0], P.VRT[2], P.VRT[3]));
    Triangle T3 = Triangle(vec3[3](P.VRT[0], P.VRT[3], P.VRT[4]));
    Triangle T4 = Triangle(vec3[3](P.VRT[0], P.VRT[4], P.VRT[1]));
    Triangle T5 = Triangle(vec3[3](P.VRT[1], P.VRT[2], P.VRT[3]));
    Triangle T6 = Triangle(vec3[3](P.VRT[3], P.VRT[4], P.VRT[1]));
    
    Triangle TR[6] = Triangle[6](T1, T2, T3, T4, T5, T6);
    vec3 tN;
    float tri;
    float t = INF;
    
    for (int i = 0; i < 6; i++) 
    {
        tri =TraceTriangle(pos, dir, TR[i], tN);
        if (tri < t) 
        {
            t = tri;
            N = tN;
        }
    }
    return t;
}

float TraceSphere(vec3 camPos, vec3 dir, Sphere S, out vec3 N) 
{
    vec3 pos = camPos - S.POS;
    float a = dot(dir, dir);
    float b = dot(pos, dir);
    float c = dot(pos, pos) - S.RAD * S.RAD;
    float D = b * b - a * c;
    if (D < 0.0) {
        return INF;
    }
    float t = -b - sqrt(D);
    if (t > 0.0) {
        N = normalize(pos + t * dir);
        return t;
    }
    t = -b + sqrt(D);
    if (t < 0.0) {
        return INF;
    }
    N = normalize(pos + t * dir);
    return t;
}

float TracePod(vec3 pos, vec3 dir, out vec3 N)
{
    float t = (-1.0 - pos.y) / dir.y;
    const float r = 2.0;
    if(t <= EPS) return INF;
    vec3 worldPos = t * dir + pos;
    if(dot(worldPos.xz, worldPos.xz) < r*r)
    {
        N = vec3(0, 1, 0);
        return t;
    }    
    float a = dot(dir.xz, dir.xz);
    float b = dot(pos.xz, dir.xz);
    float c = dot(pos.xz, pos.xz) - r*r;
    float D = b * b - a * c;    
    if(D < 0.0) return INF;
    t = (-b - sqrt(D)) / a;   
    if(t > 0.0)
    {
        worldPos = t * dir + pos;
        if (worldPos.y <= -1.0) 
        {
            N = normalize(vec3(worldPos.x, 0, worldPos.z));
            return t;
        }
    }    
    t = (-b + sqrt(D)) / a;
    if(t < EPS) return INF;
    worldPos = t * dir + pos;
    if(worldPos.y <= -1.0)
    {
        N = normalize(vec3(worldPos.x, 0, worldPos.z));
        return t;
    }
    return INF;
}

float TraceScene(vec3 pos, vec3 dir, out vec3 N, out MatType M) //  
{
    vec3 tN;
    float t = INF;
    
    float sT = TraceSphere(pos, dir, SC.S, tN);
    if (sT < t) 
    {
        t = sT;
        N = tN;
        M = SC.S.M;
    }
    
    float pT = TracePyramid(pos, dir, SC.P, tN);
    if (pT < t) 
    {
        t = pT;
        N = tN;
        M = SC.P.M;
    }
    
    vec3 worldPos = t * dir + pos;
    float cT = TracePod(pos, dir, tN);
    MatType Diff = MatType(0.0, 1.0, 0.0, 0.0, AIR_N, vec3(0));
    if(cT < t)
    {
        t = cT;
        N = tN;
        worldPos = t * dir + pos;
        M = Diff;
        M.COL = texture(iChannel3, worldPos.xz * 0.1).rgb;
    }
    
    return t;
}

bool IsOccluded(vec3 pos, LightSrc L) 
{
    vec3 dir = L.SPH.POS - pos + randDir * L.SPH.RAD;
    float dist = length(dir);
    dir /= dist;
    vec3 tN;
    MatType tM;
    float t = TraceScene(pos + dir * 1e-3, dir, tN, tM);
    return t < dist;
}

vec3 ComputeLight(vec3 pos, vec3 dir, vec3 N, MatType M) 
{
    vec3 Power = vec3(0);
    for (int i = 0; i < 3; i++) 
    {
        vec3 Light = SC.L[i].SPH.POS - pos;
        float occ = IsOccluded(pos, SC.L[i]) ? 0.5 : SC.L[i].POW / dot(Light, Light);
        Power += max(0.0, dot(N, normalize(Light))) * occ * SC.L[i].SPH.M.COL;
    }
    Power += texture(iChannel1, N).rgb * 0.1;
    return M.COL * Power;
}

vec3 Refraction(vec3 v, vec3 N, float n1, float n2) 
{
    if (dot(v, N) < 0.0) 
    {
        N = -N;
    }
    float cosA = dot(v, N);
    float sinA = sqrt(1.0 - cosA * cosA);
    vec3 tang = normalize(v - cosA * N);
    float sinB = sinA / n2 * n1;
    if (sinB > 1.0) 
    {
        return reflect(v, N);   
    }
    float cosB = sqrt(1.0 - sinB * sinB);
    return sinB * tang + cosB * N;
}

vec3 MainRender(vec3 curPos, vec3 curDir, float n1, float n2, LightSrc LS[3], vec3 randVec) 
{
    vec3 col = vec3(0.0);
    for (int j = 0; j < 8; j++) 
    {
        float t = INF;
        bool INT;
        vec3 N;
        vec3 tN;
        MatType MAT;
        
        for (int i = 0; i < 3; i++) 
        {
            float lT = TraceSphere(curPos, curDir, LS[i].SPH, tN);
            if (lT < t) 
            {
                t = lT;
                N = N;
                MAT = LS[i].SPH.M;
            }
        }
        
        MatType tM;
        float scT = TraceScene(curPos, curDir, tN, tM);
        if (scT < t) 
        {
            t = scT;
            N = tN;
            MAT = tM;
        }
        
        vec3 worldPos = t * curDir + curPos;
        float pT = TracePlane(curPos, curDir, tN);
        if (pT < t) 
        {
            t = pT;
            N = tN;
            worldPos = t * curDir + curPos;
            MAT = Diffusive;
            MAT.COL = texture(iChannel0, worldPos.xz * 0.1).rgb;
            if(randVec.y < 0.3){
                MAT = Reflective;
            }
        }
        
        if (MAT.EMIS != 0.0) 
        {
            col +=  MAT.COL * MAT.EMIS * K;
        }
        else if (MAT.DIFF != 0.0) 
        {
            col += ComputeLight(worldPos, curDir, N, MAT) * MAT.DIFF * K;
        }
        
        
        if (MAT.REFL > 0.0) {
            n1 = n2;
            n2 = MAT.N;
        }
       
        float R = powf((n1-n2)/(n1+n2), 2);
        if ((MAT.REFL != 0.0 && MAT.REFR != 0.0 && randVec.x * 0.2 < R) ||
            ((MAT.REFL == 0.0 || MAT.REFR == 0.0) && MAT.REFL != 0.0)) 
        {
            curDir = reflect(curDir,N);
            curPos = worldPos + curDir * 1e-4;
        }
        else if (( (MAT.REFL != 0.0 && MAT.REFR != 0.0) && !((MAT.REFL != 0.0 && 
            MAT.REFR != 0.0 && randVec.x * 0.2 < R) || ((MAT.REFL == 0.0 || 
            MAT.REFR == 0.0) && MAT.REFL != 0.0)) ) || ((MAT.REFL == 0.0 || 
            MAT.REFR == 0.0) &&  MAT.REFR != 0.0) )
        {
            curDir = Refraction(curDir, N, n1, n2);
            curPos = worldPos + curDir * 1e-4;
            
            if (n1 == AIR_N && n2 == GLASS_N) 
            {
                vec4 p = FireRayMarch(curPos + vec3(0.0, 2.8, 0.0), curDir);
                float glow = p.w;
                vec3 flamecol = mix( vec3(0.1, 0.5, 1), vec3(1, 0.5, 0.1), p.y * 0.022 + 2.5) ;
                col += mix(vec3(0), flamecol, pow(glow * 2.0, 4.0));
            }
        }
            
        if (t == INF) 
        {
            col += K * texture(iChannel1, curDir).rgb;
        }
    }
    return col;
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

void mainImage(out vec4 fragColor, in vec2 fragCoord )
{
    vec3 randVec = vec3(rand(float(iFrame)), rand(float(iFrame +5)), rand(float(iFrame + 15)));
    
    randDir = normalize(randVec - 0.5);
    
    MatType L1_M = Emissive;
    L1_M.COL = LIGHT_COL_1;
    MatType L2_M = Emissive;
    L2_M.COL = LIGHT_COL_2;
    MatType L3_M = Emissive;
    L3_M.COL = LIGHT_COL_3;
    
    vec3 BLUR_LIGHT_POS_1 = LIGHT_POS_1 + vec3(1,0,0)*randVec.x*0.1 
        + vec3(0,1,0)*randVec.y*0.1 + vec3(0,0,1)*randVec.z*0.1;
    
    LightSrc L1 = LightSrc(Sphere(BLUR_LIGHT_POS_1, LIGHT_RAD_1, L1_M), LIGHT_POW_1);
    LightSrc L2 = LightSrc(Sphere(LIGHT_POS_2, LIGHT_RAD_2, L2_M), LIGHT_POW_2);
    LightSrc L3 = LightSrc(Sphere(LIGHT_POS_3, LIGHT_RAD_3, L3_M), LIGHT_POW_3);
    
    LightSrc LS[3] = LightSrc[3](L1, L2, L3);
    
    Sphere S = Sphere(vec3(3, -0.2, 2.5), 1.0, Reflective);
    
    float ViewAng = 0.4;
    mat2 rot = mat2(cos(ViewAng), sin(ViewAng), -sin(ViewAng), cos(ViewAng));
    vec2 A = rot * vec2(1.5, -1.5);
    vec2 B = rot * vec2(1.5, 1.5);
    vec2 C = rot * vec2(-1.5, 1.5);
    vec2 D = rot * vec2(-1.5, -1.5);
    
    Pyramid P = Pyramid(vec3[5](vec3(0., 2., 0.), vec3(A.x, -1.99, A.y),
            vec3(B.x, -1.99, B.y), vec3(C.x, -1.99, C.y),
            vec3(D.x, -1.99, D.y)), MyGlass);
            
    SC.L = LS;
    SC.S = S;
    SC.P = P;
    
    vec2 uv = (fragCoord - iResolution.xy * 0.5 + (randVec.xy - 0.5) * 2.0)/iResolution.x;
    vec3 front = normalize(vec3(-CAMERA_POS));
    vec3 up = vec3(0,1,0);
    vec3 right = normalize(cross(front, up));
    up = normalize(cross(right, front));
    vec3 viewVec = normalize(front + right * uv.x + up* uv.y);
    vec3 col = vec3(0);
    
    col = MainRender(CAMERA_POS, viewVec, AIR_N, AIR_N, LS, randVec);
        
    fragColor = vec4(col, 1.0);
             
    
}
