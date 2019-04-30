#include <GL/gl.h>
#include <GL/glu.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <GL/glut.h>   // The GL Utility Toolkit (Glut) Header
#include <iostream>
#include <utility>

static constexpr auto WIDTH = 500;
static constexpr auto HEIGHT =  500;
static constexpr auto TO_BE_SET = 99999;
/*
vector in 3d space
*/

template <typename T>
struct vec3
{
	T x;
	T y;
	T z;
	vec3() : x(0.0), y(0.0), z(0.0) {}
	vec3(T _x, T _y, T _z) : x(_x), y(_y), z(_z) {}
	vec3(T i) : x(i), y(i), z(i) {}

	inline T dot(const vec3<T>& v) const
	{
		return x * v.x + y * v.y + z * v.z;
	}

	inline vec3<T>& normalize() noexcept
	{
		T n = x * x + y * y + z * z;
		if (n > 0)
		{
			T invN = 1 / sqrt(n);
			x *= invN;
			y *= invN;
			z *= invN;
		}
		return *this;
	}

	vec3& operator = (const vec3& rhs)
	{
		this->x = rhs.x;
		this->y = rhs.y;
		this->z = rhs.z;
		return *this;
	}

	vec3<T> operator * (const T& rhs) const
	{
		return vec3<T>  (
						x * rhs,
						y * rhs,
						z * rhs
						);
	}

	vec3<T> operator * (const vec3<T>& rhs) const
	{
		return vec3<T>  (
						x * rhs.x,
						y * rhs.y,
						z * rhs.z
						);
	}

	vec3<T> operator + (const vec3<T>& rhs) const
	{
		return vec3<T>  (
						x + rhs.x,
						y + rhs.y,
						z + rhs.z
						);
	}

	vec3<T>& operator += (const vec3& rhs)
	{
		this->x += rhs.x;
		this->y += rhs.y;
		this->z += rhs.z;
		return *this;
	}

	vec3<T> operator - (const vec3& rhs) const
	{
		return vec3<T>  (
						x - rhs.x,
						y - rhs.y,
						z - rhs.z
						);
	}

	vec3<T> operator - () const
	{
		return vec3<T>(-x, -y, -z);
	}
};

/*
Representation of a point
*/
struct Point
{
	int x;
	int y;
	vec3<float> c; //color
	Point() : x(1), y(1) {}
	Point(int _x, int _y) : x(_x), y(_y) {}
	Point(int _x, int _y, vec3<float> _c) : x(_x), y(_y), c(_c.x, _c.y, _c.z) {}
	Point& operator=(const Point& rhs)
	{
		this->x = rhs.x;
		this->y = rhs.y;
		return *this;
	}
	bool operator!=(const Point& rhs)
	{
		if (this->x != rhs.x || this->y != rhs.y)
			return true;
		return false;
	}
	bool operator<(const Point& rhs) const
	{
		if (this->x < rhs.x) return true;
		else if (this->x > rhs.x) return false;
		else return (this->y < rhs.y);
	}
	bool operator>(const Point& rhs) const
	{
		if (this->x > rhs.x) return true;
		else if (this->x < rhs.x) return false;
		else return (this->y > rhs.y);
	}
	Point operator * (const Point& rhs)
	{
		this->x *= rhs.x;
		this->y *= rhs.y;
		return *this;
	}
	inline void replacePoints(const int& _x, const int& _y) noexcept
	{
		this->x = _x;
		this->y = _y;
	}
};


void init_window()
                 /* Clear the image area, and set up the coordinate system */
{

        					       /* Clear the window */
        glClearColor(0.0,0.0,0.0,0.0);
		glShadeModel(GL_SMOOTH);
        glOrtho(0,WIDTH,0,HEIGHT,-1.0,1.0);
}

/***************************************************************************/

void write_pixel(int x, int y, vec3<float> c)
                                         /* Turn on the pixel found at x,y */
{

        //glColor3f (intensity, intensity, intensity);    
		glColor3f(c.x, c.y, c.z);             
        glBegin(GL_POINTS);
           glVertex3i( x, y, 0);
        glEnd();	
}

struct phongProperties
{
	phongProperties() : diffuseColor(0.2) {}
	//default
	float kd = .3;
	float ks = .2;
	vec3<float> diffuseColor;
	float specularExponent = 3;
};

class Circle final
{
	/*Directional Components*/
	float radius;
	float radius2;
    
public:

	/*Color Components*/
	bool transparency = false;
	bool reflect = false;

	/*Light*/
	phongProperties phong;

	vec3<float> center;
	
	/*Color Components*/
	vec3<float> surfaceColor; //color of object
	vec3<float> emissionColor; //what color light it emits

	constexpr Circle() = delete;
	~Circle() = default;
	
	/*use R value constructor*//*
	explicit Circle(
					 const vec3<float>&& _center,
					 const float&&       _radius,
					 const vec3<float>&& _surfaceColor, //color of object
					 const vec3<float>&& _emissionColor, //what color light it emits
					 const float&&       _transparency = 0.0,
					 const float&&       _reflection = 0.0
				   ) noexcept
				   : center       (std::exchange(center, _center)),
				     radius       (std::exchange(radius, _radius)),
					 radius2	  (_radius * _radius),
					 surfaceColor (std::exchange(surfaceColor, _surfaceColor)),
					 emissionColor(std::exchange(emissionColor, _emissionColor)),
					 transparency (std::exchange(transparency, _transparency)),
					 reflection   (std::exchange(reflection, _reflection))
	{}
	*/
	explicit Circle(
					 const vec3<float>& _center,
					 const float&       _radius,
					 const vec3<float>& _surfaceColor, //color of object
					 const vec3<float>& _emissionColor, //what color light it emits
					 const bool&        _transparency = false,
					 const bool&		_reflect = false
				   ) 
				   : center       (_center),
				     radius       (_radius),
					 radius2	  (_radius * _radius),
					 surfaceColor (_surfaceColor),
					 emissionColor(_emissionColor),
					 transparency (_transparency),
					 reflect	  (_reflect)
	{}

	inline constexpr float radiusSquared() noexcept { return radius*radius; }
	
	/*Geometric solution to ray-sphere intersection*/
	inline bool isIntersect(const vec3<float>& rayOrigin, const vec3<float>& rayDirection, float& t0, float t1)
	{
		vec3<float> l = center - rayOrigin;
		float tca = l.dot(rayDirection);
		if (tca < 0) return false;
		float d2 = l.dot(l) - tca * tca;
		//std::cout << d2 << " " << radius2 << "\n";
		if (d2 > radius2) return false;
		float thc = sqrt(radius2 - d2);
		t0 = tca - thc;
		t1 = tca + thc;
		
		return true; 
	}
};

struct Light
{
	explicit Light(const vec3<float>& pos, const float& inten) noexcept : position(pos), intensity(inten) {}
	vec3<float> position;
	float intensity;
};

class RayTracer final
{
	struct Objects
	{
		std::vector<Circle> circles;
		std::vector<Light> lights;

		//add objects here
		explicit Objects() noexcept 
		: circles() 
		{
			circles.emplace_back(std::move(Circle(vec3<float>(0, 0, -20), 3, vec3<float>(.5,.2,.1), vec3<float>(.7, .2, .8), false, false)));
			circles.emplace_back(std::move(Circle(vec3<float>(2, 8, -30), 4, vec3<float>(.7,.4,.5), vec3<float>(.4, .1, .8), true, false)));
			circles.emplace_back(std::move(Circle(vec3<float>(2, -8, -30), 4, vec3<float>(.2,.4,.3), vec3<float>(.1, .1, .2), false, true)));
			lights.emplace_back(std::move(Light(vec3<float>(-20, 70, -50), .7)));
		}
	};

	/*
	The Trace Function
	*/
	static bool trace(const vec3<float>& raypos, const vec3<float>& raydir, int& objectHit, float& tnear)
	{
		bool isHit = false;

		//for each object in scene
		for(auto i = 0; i < objects.circles.size(); ++i)
		{
			float t0 = TO_BE_SET;
			float t1 = TO_BE_SET;
		
			if (objects.circles[i].isIntersect(raypos, raydir, t0, t1))
			{
				if (t0 < 0) t0 = t1;
				if (t0 < tnear)
				{
					tnear = t0; //new closest object
					objectHit = i;
					isHit = true;
				}
			}
		}
		return isHit;
	}

	static vec3<float> castRay(const vec3<float>& raypos, const vec3<float>& raydir, int depth)
	{
		//check recursion depth
		if (depth > 8) return vec3<float>(.2, .2, .2);

		int objectHit = 0; //which object is hit if one is ..
		float tnear = TO_BE_SET;

		if (!trace(raypos, raydir, objectHit, tnear)) return vec3<float>(.2, .2, .2); //not hit ... return background
		
		
		/*init hit vars*/
		vec3<float> surfaceColor = objects.circles[objectHit].surfaceColor; //color of object
		vec3<float> pointHit = raypos + raydir * tnear;
		vec3<float> normalHit = pointHit - objects.circles[objectHit].center;
		normalHit.normalize(); //normalize direction
		vec3<float> hitColor = (0);


		/*init light vars*/
		vec3<float> lightAmt(0);
		vec3<float> specularColor(0);
		vec3<float> shadowOrigin = (raydir.dot(normalHit) < 0) ? pointHit + normalHit : pointHit - normalHit;
		
		/*Shade twice bc it looks good*/
		//loop over all lights
		for (auto i = 0; i < objects.lights.size(); ++i)
		{
			//init
			float tnearShadow = TO_BE_SET;
			int objectHitShadow = 0;

			//get direction -> cast light ray
			vec3<float> lightDir = objects.lights[i].position - pointHit;
			float lightDistSquared = lightDir.dot(lightDir);
			lightDir = lightDir.normalize();
			float lightNormal = std::max(0.0f, lightDir.dot(normalHit));

			//reflect light
			vec3<float> reflectionDir = lightDir - (normalHit * (2 * lightDir.dot(normalHit)));

			//calc color of specular
			specularColor += powf(std::max(0.0f, -(reflectionDir.dot(raydir))), objects.circles[objectHit].phong.specularExponent) * objects.lights[i].intensity;
		}
		
		//phong equation
		hitColor = lightAmt * objects.circles[objectHit].phong.kd * objects.circles[objectHit].phong.diffuseColor + specularColor * objects.circles[objectHit].phong.ks;
		//add color back
		hitColor = hitColor + objects.circles[objectHit].surfaceColor;

		/*Reflect*/
		if (objects.circles[objectHit].reflect)
		{
			hitColor = objects.circles[objectHit].surfaceColor;
			//reflect light
			vec3<float> reflectionDir = raydir - (normalHit * (2 * raydir.dot(normalHit)));
			//recurse, only reflect 70% of light
			hitColor += castRay(pointHit + normalHit, reflectionDir, depth+1) * 0.9f; 
		}

		/*See thru*/
		else if (objects.circles[objectHit].transparency)
		{
			//set to background to make it transparent :)
			hitColor = vec3<float>(.2, .2, .2);
			//reflect light
			vec3<float> reflectionDir = raydir - (normalHit * (2 * raydir.dot(normalHit)));
			//recurse, only reflect 70% of light
			hitColor += castRay(pointHit + normalHit, reflectionDir, depth+1) * 0.8f; 
			//hitColor += objects.circles[objectHit].surfaceColor;
		}

		/*Default phong shade*/
		else
		{
			//loop over all lights
			for (auto i = 0; i < objects.lights.size(); ++i)
			{
				//init
				float tnearShadow = TO_BE_SET;
				int objectHitShadow = 0;

				//get direction -> cast light ray
				vec3<float> lightDir = objects.lights[i].position - pointHit;
				float lightDistSquared = lightDir.dot(lightDir);
				lightDir = lightDir.normalize();
				float lightNormal = std::max(0.0f, lightDir.dot(normalHit));

				//reflect light
				vec3<float> reflectionDir = lightDir - (normalHit * (2 * lightDir.dot(normalHit)));

				//calc color of specular
				specularColor += powf(std::max(0.0f, -(reflectionDir.dot(raydir))), objects.circles[objectHit].phong.specularExponent) * objects.lights[i].intensity;
			}
			
			//phong equation
			hitColor = lightAmt * objects.circles[objectHit].phong.kd * objects.circles[objectHit].phong.diffuseColor + specularColor * objects.circles[objectHit].phong.ks;
			//add color back
			hitColor = hitColor + objects.circles[objectHit].surfaceColor;
		}

		return hitColor;
	}

public:

	//list of objects
	static Objects objects;

	/*
	The Render Function
	*/
    inline static void render()
    {
		vec3<float> pixelColor;

		float invWidth = 1 / float(WIDTH); float invHeight = 1 / float(HEIGHT);
		float fov = 30;
		float aspectratio = WIDTH / float(HEIGHT);
		float angle = tan(M_PI * 0.5 * fov / 180.);

		for(int y = 0; y < HEIGHT; ++y)
			for(int x = 0; x < WIDTH; ++x)
			{
				//do some stuff
				float xx = (2 * ((x + 0.5) * invWidth) - 1) * angle * aspectratio;
				float yy = (1 - 2 * ((y + 0.5) * invHeight)) * angle;
				vec3<float> raydir(xx, yy, -1);
				raydir.normalize(); 
				int index = 0;
				pixelColor = castRay(vec3<float>(0), raydir, index);
				write_pixel(x, y, pixelColor);
			}
    }
};

//call default constructor to make objects
RayTracer::Objects RayTracer::objects;

void display ( void )   // Create The Display Function
{

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	      // Clear Screen 
  

  // CALL YOUR CODE HERE
  RayTracer::render();

  glutSwapBuffers();                                      // Draw Frame Buffer 
}

int main(int argc, char *argv[])
{

	glutInit            ( &argc, argv ); 
    glutInitDisplayMode ( GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH ); 
	glutInitWindowSize  ( HEIGHT,WIDTH ); 
	glutCreateWindow    ( "Computer Graphics" ); 
	glutDisplayFunc     ( display );  
	glutIdleFunc	    ( display );
	
        					      
        init_window();				             //create_window
						       		
	glutMainLoop        ( ); 
    //rayTracer.render();

    return 0;
}