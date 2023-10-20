#include <stdio.h>
#define PI	3.14159265359
#define PROTON_CHARGE	1.6021892e-19
#define ELECTRIC_PERMIT_VACUUM	8.85419e-12
#define GRAVITATIONAL_CONST	6.6730e-11

typedef struct {
	double x;
	double y;
	double z;
} VECTOR;

typedef struct {
	double charge;
	double mass;
	VECTOR position;
	VECTOR force;
	VECTOR speed;
	VECTOR accel;
	VECTOR kin_E;
} OBJECT;

typedef struct {
	long p;	//amount of particles
	OBJECT* ptr; //pointer to array of particles
} SPACE;

//unit vecs
VECTOR origin = {.x=0,.y=0,.z=0};
VECTOR unit_vec_x = {.x=1, .y=0, .z=0};
VECTOR unit_vec_y = {.x=0, .y=1, .z=0};
VECTOR unit_vec_z = {.x=0, .y=0, .z=1};

//functions
//math
double rev_sqrt(double num, int iter);
//vector math
double veclenght(VECTOR a);
int veccmp(VECTOR a, VECTOR b);
VECTOR vecadd(VECTOR a, VECTOR b);
VECTOR vecsubtract(VECTOR a, VECTOR b);
VECTOR scalarprod(double scalar, VECTOR a);
VECTOR dotprod(VECTOR a, VECTOR b);
VECTOR crossprod(VECTOR a, VECTOR b);
//physics
VECTOR elec(VECTOR observation, OBJECT* particle);
VECTOR grav(VECTOR observation, OBJECT* particle);
VECTOR elecfield(VECTOR observation, SPACE* s);
VECTOR gravfield(VECTOR observation, SPACE* s);
VECTOR F_elec(VECTOR E, OBJECT* particle);
VECTOR F_grav(VECTOR G, OBJECT* particle);

//math
double rev_sqrt(double num,int iter)	{
	long i;
        double x2, y, z;
        x2 = num * 0.5F;
        y = num;
        i=*(long*)&y;
        i=0x5fe6eb50c7b537a9 - (i>>1);
        y = *(double *) &i;
        for     (i=0;i<=iter;i++)	{
                y = y * (1.5F - (x2*y*y));
        } if (y<0)	y*=-1;
        return y;
}//fast inverse square root algorithm

//vector math
double veclenght(VECTOR a)	{
	return 1/rev_sqrt(a.x*a.x+a.y*a.y+a.z*a.z,4);
}
int veccmp(VECTOR a, VECTOR b)	{
	return (a.x==b.x && a.y==b.y && a.z==b.z) ? 1 : 0;
}
VECTOR vecadd(VECTOR a, VECTOR b)	{
	return (VECTOR)	{.x=a.x+b.x,
			.y=a.y+b.y,
			.z=a.z+b.z	};
}
VECTOR vecsubtract(VECTOR a, VECTOR b)	{
	return (VECTOR)	{.x=a.x-b.x,
			.y=a.y-b.y,
			.z=a.z-b.y	};
}//note the resulting vector points from the tip of b to the tip of a
VECTOR scalarprod(double scalar, VECTOR a)	{
	return	(VECTOR) { .x=scalar*a.x,
			.y=scalar*a.y,
			.z=scalar*a.z	};
}
VECTOR dotproduct(VECTOR a, VECTOR b)	{
	return (VECTOR)	{ .x=a.x*b.x,
			.y=a.y*b.y,
			.z=a.z*b.z	};	
}//vector a is transposed
VECTOR crossproduct(VECTOR a, VECTOR b)	{
	return (VECTOR)	{ .x=(a.y*b.z-a.z*b.y),
			.y=(a.z*b.x-a.x*b.z),
			.z=(a.x*b.y-a.y*b.x)	};
}
//physics
VECTOR elec(VECTOR observation, OBJECT* particle)	{
	//E = (1/4*pi*e)*(q/r^2)
	//we measure the electrical field as result of partical at position
	//E is pointed away from the particle if q>0
	//and towards the particle if q<0
	//the direction vector is pointed away from the particle
	//so that when it is multiplied with E, and thus q>0, it points away from the particle
	VECTOR direction = vecsubtract(observation, particle->position);
	//position - particleposition gives the vector from particle
	double r = veclenght(direction);
	double E = 1/(4*PI*ELECTRIC_PERMIT_VACUUM) * (particle->charge/(r*r));
	//we normalize the direction vector, then multiply it by E
	direction = scalarprod(1/r,direction);
	direction = scalarprod(E,direction);
	return direction;
}
VECTOR grav(VECTOR observation, OBJECT* particle)	{
	//Gf=Gm/r^2
	//pointed towards particle
	VECTOR direction = vecsubtract(particle->position, observation);
	double r = veclenght(direction);
	double g = GRAVITATIONAL_CONST*particle->mass/(r*r);
	direction = scalarprod(1/r,direction);
	direction = scalarprod(g,direction);
	return direction;
}
VECTOR elecfield(VECTOR observation, SPACE* s)	{
	//elec field as result of multiple particles
	VECTOR temp=origin;
	for(int i=0;i<=(s->p);i++)	{
		if ( !veccmp( (s->ptr+i)->position,observation) )
			temp = vecadd(temp,elec(observation,(s->ptr+i)) );
	}	return temp;
}
VECTOR gravfield(VECTOR observation, SPACE* s)   {
        //elec field as result of multiple particles
        VECTOR temp=origin;
        for(int i=0;i<=(s->p);i++)        {
                if ( !veccmp( (s->ptr+i)->position,observation) )
                        temp = vecadd(temp,grav(observation,(s->ptr+i)) );
        }       return temp;
}
VECTOR F_elec(VECTOR E, OBJECT* particle)	{
	//force enacted on particle as result of electrical field
	return scalarprod(particle->charge,E);
}
VECTOR F_grav(VECTOR G, OBJECT* particle)	{
	//force enacted on particle as result of electrical field
	return scalarprod(particle->mass,G);
}
