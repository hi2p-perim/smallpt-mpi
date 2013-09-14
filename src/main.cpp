#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <mpi.h>
#include <queue>

struct Vec
{
	double x, y, z;
	Vec(double x_ = 0, double y_ = 0, double z_ = 0) { x=x_; y=y_; z=z_; }
	Vec operator+(const Vec &b) const { return Vec(x + b.x, y + b.y, z + b.z); }
	Vec operator-(const Vec &b) const { return Vec(x - b.x, y - b.y, z - b.z); }
	Vec operator*(double b) const { return Vec(x * b, y * b, z * b); }
	Vec mult(const Vec &b) const { return Vec(x * b.x, y * b.y, z * b.z); }
	Vec& norm() { return *this = *this * (1 / sqrt(x * x + y * y + z * z)); }
	double dot(const Vec &b) const { return x * b.x + y * b.y + z * b.z; } // cross:
	Vec operator%(Vec&b) { return Vec(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x); }
};

struct Ray
{
	Vec o, d; Ray(Vec o_, Vec d_)
		: o(o_), d(d_)
	{}
};

enum Refl_t
{
	DIFF,
	SPEC,
	REFR
};

struct Sphere
{
	double rad;       // radius
	Vec p, e, c;      // position, emission, color
	Refl_t refl;      // reflection type (DIFFuse, SPECular, REFRactive)
	
	Sphere(double rad_, Vec p_, Vec e_, Vec c_, Refl_t refl_)
		: rad(rad_)
		, p(p_)
		, e(e_)
		, c(c_)
		, refl(refl_)
	{}
	
	double intersect(const Ray &r) const 
	{
		// returns distance, 0 if nohit
		Vec op = p-r.o; // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
		double t, eps=1e-4, b=op.dot(r.d), det=b*b-op.dot(op)+rad*rad;
		if (det<0) return 0; else det=sqrt(det);
		return (t=b-det)>eps ? t : ((t=b+det)>eps ? t : 0);
	}
};

Sphere spheres[] =
{
	//Scene: radius, position, emission, color, material
	Sphere(1e5, Vec( 1e5+1,40.8,81.6), Vec(),Vec(.75,.25,.25),DIFF), //Left
	Sphere(1e5, Vec(-1e5+99,40.8,81.6),Vec(),Vec(.25,.25,.75),DIFF), //Rght
	Sphere(1e5, Vec(50,40.8, 1e5),     Vec(),Vec(.75,.75,.75),DIFF), //Back
	Sphere(1e5, Vec(50,40.8,-1e5+170), Vec(),Vec(),           DIFF), //Frnt
	Sphere(1e5, Vec(50, 1e5, 81.6),    Vec(),Vec(.75,.75,.75),DIFF), //Botm
	Sphere(1e5, Vec(50,-1e5+81.6,81.6),Vec(),Vec(.75,.75,.75),DIFF), //Top
	Sphere(16.5,Vec(27,16.5,47),       Vec(),Vec(1,1,1)*.999, SPEC), //Mirr
	Sphere(16.5,Vec(73,16.5,78),       Vec(),Vec(1,1,1)*.999, REFR), //Glas
	Sphere(600, Vec(50,681.6-.27,81.6),Vec(12,12,12),  Vec(), DIFF)  //Lite
};

inline double clamp(double x)
{
	return x < 0 ? 0 : x > 1 ? 1 : x;
}

inline int toInt(double x)
{
	return int(pow(clamp(x), 1 / 2.2) * 255 + .5);
}

inline bool intersect(const Ray &r, double &t, int &id)
{
	double n = sizeof(spheres) / sizeof(Sphere);
	double d;
	double inf = t = 1e20;
	
	for (int i = int(n); i--;)
	{
		if ((d=spheres[i].intersect(r))&&d<t) { t=d;id=i; }
	}
	
	return t < inf;
}

Vec radiance(const Ray &r, int depth, unsigned short *Xi)
{
	double t;	// distance to intersection
	int id=0;	// id of intersected object

	if (!intersect(r, t, id))
	{
		// if miss, return black
		return Vec(); 
	}

	// the hit object
	const Sphere &obj = spheres[id];       

	Vec x = r.o+r.d*t;
	Vec n = (x - obj.p).norm();
	Vec nl = n.dot(r.d) < 0? n : n * -1;
	Vec f = obj.c;

	// max refl
	double p = f.x > f.y && f.x > f.z ? f.x : f.y>f.z ? f.y : f.z;

	if (++depth > 5)
	{
		if (erand48(Xi) < p)
		{
			f = f * (1 / p);
		}
		else
		{
			//R.R.
			return obj.e;
		}
	}

	if (obj.refl == DIFF)
	{
		// Ideal DIFFUSE reflection
		double r1 = 2 * M_PI * erand48(Xi);
		double r2 = erand48(Xi);
		double r2s = sqrt(r2);

		Vec w = nl;
		Vec u = ((fabs(w.x) > .1 ? Vec(0,1) : Vec(1)) % w).norm();
		Vec v = w % u;
		Vec d = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).norm();

		return obj.e + f.mult(radiance(Ray(x,d), depth, Xi));
	}
	else if (obj.refl == SPEC)
	{
		// Ideal SPECULAR reflection
		return obj.e + f.mult(radiance(Ray(x, r.d - n * 2 * n.dot(r.d)), depth, Xi));
	}
	else
	{
		// Ideal dielectric REFRACTION
		Ray reflRay(x, r.d-n*2*n.dot(r.d));

		// Ray from outside going in?
		bool into = n.dot(nl)>0;
		double nc = 1;
		double nt = 1.5;
		double nnt = into ? nc / nt : nt / nc;
		double ddn = r.d.dot(nl);
		double cos2t;

		if ((cos2t = 1 - nnt * nnt * (1 - ddn * ddn)) < 0)
		{
			// Total internal reflection
			return obj.e + f.mult(radiance(reflRay, depth, Xi));
		}

		Vec tdir = (r.d * nnt - n * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t)))).norm();
		double a = nt - nc;
		double b = nt + nc;
		double R0 = a * a / (b * b);
		double c = 1 - (into ? -ddn : tdir.dot(n));
		double Re = R0 + (1 - R0) * c * c * c * c * c;
		double Tr = 1 - Re;
		double P = .25 + .5 * Re;
		double RP = Re / P;
		double TP = Tr / (1 - P);

		return
			obj.e + f.mult(
				// Russian roulette
				depth > 2
					? (erand48(Xi) < P ? radiance(reflRay,depth,Xi) * RP : radiance(Ray(x, tdir), depth, Xi) * TP)
					: radiance(reflRay, depth, Xi) * Re + radiance(Ray(x, tdir), depth, Xi) * Tr);
	}
}

enum TagType
{
	TagType_Result = 0,
	TagType_AssignTask = 1,
	TagType_Exit = 2
};

const int w = 1024;
const int h = 768;

struct TaskResult
{
	int y;			// Column number
	Vec data[w];	// Data
};

int main(int argc, char *argv[])
{
	int rank;
	int num_procs;
	int proc_name_len;
	char proc_name[MPI_MAX_PROCESSOR_NAME];

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Get_processor_name(proc_name, &proc_name_len);

	// Define struct
	int blocklens[2] = { 1, w * 3 };
	MPI_Aint indices[2] = { 0, 4 };
	MPI_Datatype type[2] = { MPI_INT, MPI_DOUBLE };
	MPI_Datatype result_type;
	MPI_Type_struct(2, blocklens, indices, type, &result_type);
	MPI_Type_commit(&result_type);

	int samps = argc == 2 ? atoi(argv[1]) / 4 : 1; // # samples
	Vec* c = new Vec[w * h];
	TaskResult* result = new TaskResult;
	MPI_Status status;

	if (rank == 0)
	{
		// Master process

		// First assign initial tasks to slave processes
		for (int i = 1; i < num_procs; i++)
		{
			int y = i - 1;
			MPI_Send(&y, 1, MPI_INT, i, TagType_AssignTask, MPI_COMM_WORLD);
		}

		// Each task executes rendering of single column
		for (int y = num_procs - 1; y < h; y++)
		{
			// Wait for result
			MPI_Recv(result, 1, result_type, MPI_ANY_SOURCE, TagType_Result, MPI_COMM_WORLD, &status);
			memcpy(&c[result->y * w].x, &result->data, sizeof(double) * 3 * w);

			// Now send assigned task id (column index) to the slave process
			MPI_Send(&y, 1, MPI_INT, status.MPI_SOURCE, TagType_AssignTask, MPI_COMM_WORLD);
		}

		// There is no more tasks
		// Receive the remaining results
		for (int i = 1; i < num_procs; i++)
		{
			MPI_Recv(result, 1, result_type, MPI_ANY_SOURCE, TagType_Result, MPI_COMM_WORLD, &status);
			memcpy(&c[result->y * w].x, &result->data, sizeof(double) * 3 * w);
		}

		// Send message to exit
		for (int i = 1; i < num_procs; i++)
		{
			MPI_Send(NULL, 0, MPI_INT, i, TagType_Exit, MPI_COMM_WORLD);
		}

		// Write image to PPM file.
		FILE *f = fopen("image.ppm", "w");
		fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
		for (int i = 0; i < w * h; i++)
		{
			fprintf(f, "%d %d %d ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z));
		}
	}
	else
	{
		// Slave process
		int y;

		while (true)
		{
			// Receive task assignment
			MPI_Recv(&y, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	
			if (status.MPI_TAG == TagType_Exit)
			{
				std::cout << "Finished " << proc_name << std::endl;
				break;
			}

			//std::cout << "Processing y = " << y << " in " << proc_name << std::endl;

			Ray cam(Vec(50, 52, 295.6), Vec(0, -0.042612, -1).norm()); // cam pos, dir
			Vec cx = Vec(w * .5135 / h);
			Vec cy = (cx % cam.d).norm() * .5135;
			Vec r;

			unsigned short Xi[3] = { 0, 0, y * y * y };

			#pragma omp parallel for schedule(dynamic, 1) private(r)
			for (int x = 0; x < w; x++)
			{
				result->data[x] = Vec();

				// Loop cols
				for (int sy = 0; sy < 2; sy++)
				{
					// 2x2 subpixel rows
					for (int sx = 0; sx < 2; sx++, r = Vec())
					{
						// 2x2 subpixel cols
						for (int s = 0; s < samps; s++)
						{
							double r1 = 2 * erand48(Xi);
							double dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
							double r2 = 2 * erand48(Xi);
							double dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);

							Vec d =
								cx * (((sx + .5 + dx) / 2 + x) / w - .5) +
								cy * (((sy + .5 + dy) / 2 + y) / h - .5) + cam.d;

							r = r + radiance(Ray(cam.o + d * 140, d.norm()), 0, Xi) * (1. / samps);
						}

						result->data[x] = result->data[x] + Vec(clamp(r.x), clamp(r.y), clamp(r.z)) * .25;
					}
				}
			}

			// Send result
			result->y = h - y - 1;
			MPI_Send(result, 1, result_type, 0, TagType_Result, MPI_COMM_WORLD);
		}
	}

	MPI_Finalize();

	return 0;
}

