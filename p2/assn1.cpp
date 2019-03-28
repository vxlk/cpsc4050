//
//		          Programming Assignment #1 
//
//			        Victor Zordan
//		
//		
//
/***************************************************************************/

                                                   /* Include needed files */
#include <GL/gl.h>
#include <GL/glu.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <GL/glut.h>   // The GL Utility Toolkit (Glut) Header
#include <iostream>
#include <map>
#include <time.h>
#include "objParser.h"

#define WIDTH 1000
#define HEIGHT 1000

//matrix class

template<typename T>
class Matrix44
{
public:

T x[4][4] = {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};

Matrix44() {}

Matrix44 (T a, T b, T c, T d, T e, T f, T g, T h,
T i, T j, T k, T l, T m, T n, T o, T p)
{
x[0][0] = a;
x[0][1] = b;
x[0][2] = c;
x[0][3] = d;
x[1][0] = e;
x[1][1] = f;
x[1][2] = g;
x[1][3] = h;
x[2][0] = i;
x[2][1] = j;
x[2][2] = k;
x[2][3] = l;
x[3][0] = m;
x[3][1] = n;
x[3][2] = o;
x[3][3] = p;
}

const T* operator [] (int i) const { return x[i]; }
T* operator [] (int i) { return x[i]; }

Matrix44 operator * (const Matrix44& v) const
{
	Matrix44 tmp;
	for (int i = 0; i < 4; ++i) 
	{
		for (int j = 0; j < 4; ++j) 
		{
		tmp[i][j] = this->x[i][0] * v[0][j] + this->x[i][1] * v[1][j] +
				    this->x[i][2] * v[2][j] + this->x[i][3] * v[3][j];
		}
	} 

	return tmp;
} 

Matrix44 operator = (const Matrix44& v)
{
	for (int i = 0; i < 4; ++i)
		for (int j = 0; j < 4; ++j)
			this->x[i][j] = v.x[i][j];

	return *this;
}

};

int x_last,y_last;
int vertSize;

template <typename T>
struct vec3
{
	T x;
	T y;
	T z;
	vec3() : x(0.0), y(0.0), z(0.0) {}
	vec3(T _x, T _y, T _z) : x(_x), y(_y), z(_z) {}
	vec3(T i) : x(i), y(i), z(i) {}

	vec3 operator * (const vec3& rhs)
	{
		this->x *= rhs.x;
		this->y *= rhs.y;
		this->z *= rhs.z;
		return *this;
	}

	vec3 operator + (const vec3& rhs)
	{
		this->x += rhs.x;
		this->y += rhs.y;
		this->z += rhs.z;
		return *this;
	}
};

//used to store transformations
struct Transformations
{
	static void reset()
	{
		rotationMatrix = blank;
		scaleMatrix = blank;
		translationMatrix = blank;
	}
	// i think this avoids gimbal lock
	static Matrix44<double> rotationMatrix;
	// save the current scale
	static Matrix44<double> scaleMatrix;
	// save the current translation
	static Matrix44<double> translationMatrix;

	static Matrix44<double> blank;
};

Matrix44<double> scaleM;
Matrix44<double> Transformations::blank;
Matrix44<double> Transformations::rotationMatrix = blank;
Matrix44<double> Transformations::scaleMatrix = blank;
Matrix44<double> Transformations::translationMatrix = blank;

bool rotation = false;
bool translation = false;

bool scale = false;
bool scaleUp = false;
bool scaleDown = false;

bool translateW = false;
bool translateA = false;
bool translateS = false;
bool translateD = false;

bool rotateW = false;
bool rotateA = false;
bool rotateS = false;
bool rotateD = false;

std::vector<vec3<double> > verts;
std::vector<vec3<double> > vertsOriginal; //to restore points
std::vector<int> faces;

//decision variables for rendering
bool changed = true;

bool orthographic = true;
bool perspective = false;

//utility for perspective proj
void setScreenCoordsPersp(
					const double &angleOfView,
					const double &imageAspectRatio,
					const double &n, const double &f,
					double &b, double &t, double &l, double &r
				   )
	{
	float scale = tan(angleOfView * 0.5 * M_PI / 180) * n;
	r = imageAspectRatio * scale, l = -r;
	t = scale, b = -t;
} 

//perspective projection normalization
void perspProj(const float &b, const float &t, const float &l, const float &r,
				const float &n, const float &f,
				Matrix44<double> &M)
{
	M[0][0] = 2 * n / (r - l);
	M[0][1] = 0;
	M[0][2] = 0;
	M[0][3] = 0;

	M[1][0] = 0;
	M[1][1] = 2 * n / (t - b);
	M[1][2] = 0;
	M[1][3] = 0;

	M[2][0] = (r + l) / (r - l);
	M[2][1] = (t + b) / (t - b);
	M[2][2] = -(f + n) / (f - n);
	M[2][3] = -1;

	M[3][0] = 0;
	M[3][1] = 0;
	M[3][2] = -2 * f * n / (f - n);
	M[3][3] = 0; 
} 

/*
Orthographic Projection normalization
*/
void orthoProj(
				const double &b, const double &t, const double &l, 
				const double &r, const double &f, const double &n,
				Matrix44<double> &M
			  ) 
{
	M[0][0] = 2 / (r - l);
	M[0][1] = 0;
	M[0][2] = 0;
	M[0][3] = 0;

	M[1][0] = 0;
	M[1][1] = 2 / (t - b);
	M[1][2] = 0;
	M[1][3] = 0;

	M[2][0] = 0;
	M[2][1] = 0;
	M[2][2] = -2 / (f - n);
	M[2][3] = 0;

	M[3][0] = -(r + l) / (r - l);
	M[3][1] = -(t + b) / (t - b);
	M[3][2] = -(f + n) / (f - n);
	M[3][3] = 1; 
}

void multPointMatrix(const vec3<double>& in, const Matrix44<double>& M, vec3<double>& out)
{
	//out = in * Mproj;
	out.x = in.x * M[0][0] + in.y * M[1][0] + in.z * M[2][0] + /* in.z = 1 */ M[3][0];
	out.y = in.x * M[0][1] + in.y * M[1][1] + in.z * M[2][1] + /* in.z = 1 */ M[3][1];
	out.z = in.x * M[0][2] + in.y * M[1][2] + in.z * M[2][2] + /* in.z = 1 */ M[3][2];
	double w = in.x * M[0][3] + in.y * M[1][3] + in.z * M[2][3] + /* in.z = 1 */ M[3][3];

	//normalize
	if (w != 1) 
	{
	out.x /= w;
	out.y /= w;
	out.z /= w;
	} 
}

void translate(Matrix44<double>& M, const double& x = 1.0, const double& y = 1.0, const double& z = 1.0)
{
	Matrix44<double> tMatrix(1,0,0,0, 0,1,0,0, 0,0,1,0, x,y,z,1);
	//M[3][0] *= x;
	//M[3][2] *= z;
	Transformations::translationMatrix = Transformations::translationMatrix * tMatrix;
	//Transformations::translationMatrix[0][3] *= x;
	//Transformations::translationMatrix[1][3] *= y;
	//Transformations::translationMatrix[2][3] *= z;
}

void scaleMatrixF(Matrix44<double>& M, bool down = false)
{
	if (down)
	{
		Matrix44<double> downMatrix(.9,0,0,0, 0,.9,0,0, 0,0,.9,0, 0,0,0,1);
		/*
		M[0][0] /= 1.5;
		M[1][1] /= 1.5;
		M[2][2] /= 1.5;
		*/
		/*
		Transformations::scaleMatrix[0][0] /= 1.1;
		Transformations::scaleMatrix[1][1] /= 1.1;
		Transformations::scaleMatrix[2][2] /= 1.1;
		*/
		Transformations::scaleMatrix = Transformations::scaleMatrix * downMatrix;
		//scaleM = scaleM * downMatrix;
		
		
	}

	else
	{
		/*
		M[0][0] *= 1.5;
		M[1][1] *= 1.5;
		M[2][2] *= 1.5; 
		*/
		/*
		Transformations::scaleMatrix[0][0] *= 1.1;
		Transformations::scaleMatrix[1][1] *= 1.1;
		Transformations::scaleMatrix[2][2] *= 1.1;
		*/
		Matrix44<double> upMatrix(1.1,0,0,0, 0,1.1,0,0, 0,0,1.1,0, 0,0,0,1);
		Transformations::scaleMatrix = Transformations::scaleMatrix * upMatrix;
	}
	
}

void rotate(Matrix44<double>& M, bool isLeftOrDown = false, bool isHorizontalAxis = false)
{
	double amt = 1.2 *(M_PI /(double) 180.0);
	
	//rotate around horizontal
	if (isHorizontalAxis)
	{
		//rotate down -
		if (isLeftOrDown)
		{
			Matrix44<double> ld(1,0,0,0, 0,cos(-amt),-sin(-amt),0, 0,sin(-amt),cos(-amt),0, 0,0,0,1);	
			/*
			M[1][1] *= cos(-20.0);
			M[1][2] *= -sin(-20.0);
			M[2][1] *= sin(-20.0);
			M[2][2] *= cos(-20.0);
			*/
		/*
			Transformations::rotationMatrix[1][1] *= cos(-20.0);
			Transformations::rotationMatrix[1][2] *= -sin(-20.0);
			Transformations::rotationMatrix[2][1] *= sin(-20.0);
			Transformations::rotationMatrix[2][2] *= cos(-20.0);
			*/
			Transformations::rotationMatrix = Transformations::rotationMatrix * ld;
		}
		else
		{
			Matrix44<double> nld(1,0,0,0, 0,cos(amt),-sin(amt),0, 0,sin(amt),cos(amt),0, 0,0,0,1);
			/*
			Transformations::rotationMatrix[1][1] *= cos(120);
			Transformations::rotationMatrix[1][2] *= -sin(120);
			Transformations::rotationMatrix[2][1] *= sin(120);
			Transformations::rotationMatrix[2][2] *= cos(120);
			*/
			Transformations::rotationMatrix = Transformations::rotationMatrix * nld;
		}
		
	}

	//rotate around vertical
	else
	{
		//rotate left -
		if (isLeftOrDown)
		{
			Matrix44<double> dl(cos(amt),0,sin(-amt),0, 0,1,0,0, -sin(-amt),0,cos(-amt),0, 0,0,0,1);
			/*
			Transformations::rotationMatrix[0][0] *= cos(-20.0);
			Transformations::rotationMatrix[0][2] *= sin(-20.0);
			Transformations::rotationMatrix[2][0] *= -sin(-20.0);
			Transformations::rotationMatrix[2][2] *= cos(-20.0);
			*/
			Transformations::rotationMatrix = Transformations::rotationMatrix * dl;
		}

		//rotate right
		else
		{
			Matrix44<double> dln(cos(amt),0,sin(amt),0, 0,1,0,0, -sin(amt),0,cos(amt),0, 0,0,0,1);
			/*
			Transformations::rotationMatrix[0][0] *= cos(20.0);
			Transformations::rotationMatrix[0][2] *= sin(20.0);
			Transformations::rotationMatrix[2][0] *= -sin(20.0);
			Transformations::rotationMatrix[2][2] *= cos(20.0);
			*/
			Transformations::rotationMatrix = Transformations::rotationMatrix * dln;
		}
	}
	
}
/*
Representation of a color
*/
struct Color
{
	double r;
	double g;
	double b;
	Color(double _r, double _g, double _b) : r(_r), g(_g), b(_b) {}
	Color() : r(1.0), g(1.0), b(1.0) {}
	void update(double _r, double _g, double _b) { r = _r; g = _g; b = _b; }
	void randomize()
	{
		srand(time(NULL));
		this->r = (double)(rand()* (1.0 - .01)) / RAND_MAX;
		this->g = (double)(rand()* (1.0 - .01)) / RAND_MAX;
		this->b = (double)(rand()* (1.0 - .01)) / RAND_MAX;
	}
	Color operator-(const Color& rhs) const
	{
		Color toBeReturned;
		toBeReturned.r -= rhs.r;
		toBeReturned.g -= rhs.g;
		toBeReturned.b -= rhs.b;
		return toBeReturned;
	}
	Color operator+(const Color& rhs) const
	{
		Color toBeReturned;
		toBeReturned.r += rhs.r;
		toBeReturned.g += rhs.g;
		toBeReturned.b += rhs.b;
		return toBeReturned;
	}
	Color operator/(const int& rhs)
	{
		Color toBeReturned;
		toBeReturned.r = toBeReturned.r/(double)rhs;
		toBeReturned.g = toBeReturned.g/(double)rhs;
		toBeReturned.b = toBeReturned.b/(double)rhs;
		return toBeReturned;
	}
	Color operator*(const double& rhs)
	{
		Color toBeReturned;
		toBeReturned.r = toBeReturned.r*(double)rhs;
		toBeReturned.g = toBeReturned.g*(double)rhs;
		toBeReturned.b = toBeReturned.b*(double)rhs;
		return toBeReturned;
	}
	bool operator==(const Color& rhs)
	{
		return (this->r == rhs.r && this->g == rhs.g && this->b == rhs.b);
	}
	//uint8_t a;
};

//forward declare
void write_pixel(int x, int y, Color c);

/*
Representation of a point
*/
struct Point
{
	int x;
	int y;
	Color c;
	Point() : x(1), y(1) {}
	Point(int _x, int _y) : x(_x), y(_y) {}
	Point(int _x, int _y, Color _c) : x(_x), y(_y), c(_c.r, _c.g, _c.b) {}
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
	void replacePoints(const int& _x, const int& _y)
	{
		this->x = _x;
		this->y = _y;
	}
};

//a buffer of all the previous points ... now I can draw based solely on position in buffer
std::vector<Point> listOfPoints;
//to draw separate curves
std::vector<int> indexes;
//last picked color
Color selectedColor(1.0,1.0,1.0);

/*
Connects two points in space
*/
struct Line
{
	static bool isToggled;

	static void drawFace(const vec3<double>& v1, const vec3<double>& v2, const vec3<double>& v3)
	{
		draw(Point(v1.x, v1.y), Point(v2.x, v2.y));
		draw(Point(v2.x, v2.y), Point(v3.x, v3.y));
		draw(Point(v1.x, v1.y), Point(v3.x, v3.y));
	}

	static void draw(const Point& a, const Point& b)
	{
		//points to be drawn
		float x = a.x;
		float y = a.y; //b
		
		//connect old x,y -> new x,y
		//y = (Y2-Y1)/(X2-X1)x + b


		//calculate the number of steps....
		int steps = abs(b.x - a.x) > abs(b.y - a.y) ? 
					abs(b.x - a.x) : 
					abs(b.y - a.y);

		/*
		Color incrementer = b.c - a.c;
		incrementer = incrementer /(double) steps;
		*/
		Color pixelColor;
		
		//calculate increment values
		float changeX = (b.x - a.x)/(float)steps;
		float changeY = (b.y - a.y)/(float)steps;

		for (int i = 0; i <= steps; ++i)
		{
			double u = double(i)/(steps);
			pixelColor.r = a.c.r * (1.0 - u) + b.c.r * u;
			pixelColor.g = a.c.g * (1.0 - u) + b.c.g * u;
			pixelColor.b = a.c.b * (1.0 - u) + b.c.b * u;
			//std::cout << pixelColor.r << " " << pixelColor.g << " " << pixelColor.b << "\n";
			write_pixel(x,y,pixelColor); //temp intensity COLOR GOES HERE

			//decide next pixel, the casting to an int decides the rounding
			x += changeX;
			y += changeY;
		}
	}

	static void draw(const int& pos)
	{
		//points to be drawn
		float x = listOfPoints[pos].x;
		float y = listOfPoints[pos].y;
		
		//connect old x,y -> new x,y
		//y = (Y2-Y1)/(X2-X1)x + b


		//calculate the number of steps....
		int steps = abs(listOfPoints[pos-1].x - listOfPoints[pos].x) > abs(listOfPoints[pos-1].y - listOfPoints[pos].y) ? 
					abs(listOfPoints[pos-1].x - listOfPoints[pos].x) : 
					abs(listOfPoints[pos-1].y - listOfPoints[pos].y);

		
		//calculate increment values
		float changeX = (listOfPoints[pos-1].x - listOfPoints[pos].x)/(float)steps;
		float changeY = (listOfPoints[pos-1].y - listOfPoints[pos].y)/(float)steps;
		Point a = listOfPoints[pos];
		Point b = listOfPoints[pos-1];
		Color pixelColor;
		for (int i = 0; i <= steps; ++i)
		{
			double u = double(i)/(steps);
			pixelColor.r = a.c.r * (1.0 - u) + b.c.r * u;
			pixelColor.g = a.c.g * (1.0 - u) + b.c.g * u;
			pixelColor.b = a.c.b * (1.0 - u) + b.c.b * u;
			//std::cout << pixelColor.r << " " << pixelColor.g << " " << pixelColor.b << "\n";
			write_pixel(x,y,pixelColor); //temp intensity COLOR GOES HERE

			//decide next pixel, the casting to an int decides the rounding
			x += changeX;
			y += changeY;
		}
	}
};

/*
Base class for a curve
*/
struct Curve
{


};
/*
Curve using a bezier algorithm
*/
struct BezierCurve : public Curve
{
	static bool isToggled;

	static void draw(const int& pos)
	{
		//go back 4 positions to get the last 4 points
		Point p0 = listOfPoints[pos];
		Point p1 = listOfPoints[pos-1];
		Point p2 = listOfPoints[pos-2];
		Point p3 = listOfPoints[pos-3];

		Point newPoint;
		Point lastPoint = p0; //set starting point

		//decide increment

		for(float t = 0; t < 1; t+=.01)
		{
			newPoint.x = pow(1.0-t,3.0)*p0.x + 3.0*pow(1.0-t,2.0)*t*p1.x + 3.0*(1-t)*pow(t,2.0)*p2.x + pow(t,3.0)*p3.x;
			newPoint.y = pow(1.0-t,3.0)*p0.y + 3.0*pow(1.0-t,2.0)*t*p1.y + 3.0*(1-t)*pow(t,2.0)*p2.y + pow(t,3.0)*p3.y;
			Line::draw(newPoint, lastPoint);
			lastPoint = newPoint;
		}

	}
};

class ColorPicker
{
	std::map<Point, Color> palette;

public:
	static bool isToggled;
	ColorPicker() {}
	void draw()
	{
		double inc = .025;
		double r = 0.0;
		double g = 0.0;
		double b = 0.0;

		
			int x = 0;
			int y = 0;
			const int max = 250;
			while(r < 1.0)
			{
				palette.insert(std::pair<Point, Color>(Point(x,y), Color(r,g,b)));
				write_pixel(x, y, Color(r,g,b));
				/*
				if (y < max) ++y;
				else
				{
					++x;
					y = 0;
				}
				*/
				r+=inc;
				g = 0;
				while(g < 1.0)
				{
					palette.insert(std::pair<Point, Color>(Point(x,y), Color(r,g,b)));
					write_pixel(x, y, Color(r,g,b));
					/*
					if (y < max) ++y;
					else
					{
						++x;
						y = 0;
					}
					*/
					g+=inc;
					b = 0;
					while(b < 1.0)
					{
						palette.insert(std::pair<Point, Color>(Point(x,y), Color(r,g,b)));
						write_pixel(x, y, Color(r,g,b));
						if (y < max) ++y;
						else
						{
							++x;
							y = 0;
						}
						b+=inc;
					}
				}
			}
			//palette.insert(std::pair<Point, Color>(Point(i,u), Color(r,g,b)));
			//write_pixel(i, u, Color(r,g,b));
			//}

	}

	Color updateColor(Point clickedPos)
	{
		//check the point in the circle, return the color
		//selectedColor.update()
		if (this->palette.find(clickedPos)!=this->palette.end())
			return palette[clickedPos];
		else return selectedColor;
	}
};

//Color picker gui
ColorPicker cp;
bool Line::isToggled = false;
bool BezierCurve::isToggled = false;
bool ColorPicker::isToggled = false;

/***************************************************************************/

void init_window()
                 /* Clear the image area, and set up the coordinate system */
{

        					       /* Clear the window */
        glClearColor(0.0,0.0,0.0,0.0);
		glShadeModel(GL_SMOOTH);
        glOrtho(0,WIDTH,0,HEIGHT,-1.0,1.0);
}

/***************************************************************************/

void write_pixel(int x, int y, Color c)
                                         /* Turn on the pixel found at x,y */
{

        //glColor3f (intensity, intensity, intensity);    
		glColor3f(c.r, c.g, c.b);             
        glBegin(GL_POINTS);
           glVertex3i( x, y, 0);
        glEnd();	
}

//***************************************************************************/

void display ( void )   // Create The Display Function
{
	
	
   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	      // Clear Screen 

  //draw a pixel at mouse pos
 // write_pixel(x_last,y_last,selectedColor);

  // CALL YOUR CODE HERE
  if (changed)
	{
		changed = false;
		Matrix44<double> projMatrix;
		Matrix44<double> worldToCamera;
		
		//set camera to good view dist
		worldToCamera[3][1] = -20;
		worldToCamera[3][2] = -50;
		
		//set params
		double near = 1;
		double far = 90;
		double aspectRatio = WIDTH / HEIGHT;

		vec3<double> minWorld(99999.0);
		vec3<double> maxWorld(-99999.0);
		
		//compute bounding box
		
		for (int i = 0; i < vertSize; ++i)
		{
			if (verts[i].x < minWorld.x) minWorld.x = verts[i].x;
			if (verts[i].y < minWorld.y) minWorld.y = verts[i].y;
			if (verts[i].z < minWorld.z) minWorld.z = verts[i].z;
			if (verts[i].x > maxWorld.x) maxWorld.x = verts[i].x;
			if (verts[i].y > maxWorld.y) maxWorld.y = verts[i].y;
			if (verts[i].z > maxWorld.z) maxWorld.z = verts[i].z;
		}
		
		//xform min and max to camera
		vec3<double> minCamera;
		vec3<double> maxCamera;
		
		//get largest points
		double maxX = std::max(fabs(minWorld.x), fabs(maxWorld.x));
		double maxY = std::max(fabs(minWorld.y), fabs(maxWorld.y));

		double max = std::max(maxX, maxY);

		double right;
		double t;
		double l;
		double b;

		//project
		if (orthographic)
		{
			right = max * aspectRatio;
			t = max;
			l = -right;
			b = -t;
			orthoProj(b, t, l, right, near, far, projMatrix);
		}
		else if (perspective)
		{
			setScreenCoordsPersp(90.0, aspectRatio, near, far, b, t, l, right); 
			perspProj(b, t, l, right, near, far, projMatrix); 
		}
		
		// loop over all verts
		for (int i = 0; i <vertSize; ++i)
		{
			//init
			vec3<double> vertCamera;
			vec3<double> projectedVert;

			//transform
			if (translation)
			{
				if (translateW)
					translate(Transformations::translationMatrix, 1.0, 2.0);
				if (translateS)
					translate(Transformations::translationMatrix, 1.0, -2.0);
				if (translateA)
					translate(Transformations::translationMatrix, 2.0);
				if (translateD)
					translate(Transformations::translationMatrix, -2.0);
				multPointMatrix(verts[i], Transformations::translationMatrix, verts[i]);
				
			}

			if (rotation)
			{
				if (rotateW)
					rotate(projMatrix, false, true);
				if (rotateS)
					rotate(projMatrix, true, true);
				if (rotateA)
					rotate(projMatrix, true, false);
				if (rotateD)
					rotate(projMatrix, false, false);
				multPointMatrix(verts[i], Transformations::rotationMatrix, verts[i]);
			}

			if (scale)
			{
				if (scaleUp)
					scaleMatrixF(Transformations::scaleMatrix);
				else if (scaleDown)
					scaleMatrixF(Transformations::scaleMatrix, true);
				multPointMatrix(verts[i], Transformations::scaleMatrix, verts[i]);
			}

			//xform to camera space
			multPointMatrix(verts[i], worldToCamera, vertCamera);

			//project
			multPointMatrix(vertCamera, projMatrix, projectedVert);
			

			//transform back from normalized points
			verts[i].x = (projectedVert.x * WIDTH) /(int) 2;
			verts[i].y = (projectedVert.y * HEIGHT) /(int) 2;
			
		}
	}

	//reset -> keep the display from rendering an infinite tranformation
	scaleUp = false;
	scaleDown = false;
	translateW = false;
	translateA = false;
	translateS = false;
	translateD = false;
	rotateW = false;
	rotateA = false;
	rotateS = false;
	rotateD = false;

	//render
	for (int i = 0; i < faces.size(); i+=3)
	{
	//	if (!i) std::cout << "start---------------------------------------------\n";
		//faces stored as 1 index -> convert to 0 index
		Line::drawFace(verts[faces[i]-1], verts[faces[i+1]-1], verts[faces[i+2]-1]);
		std::cout << Transformations::rotationMatrix[0][0] << "\n";
		//std::cout << verts[faces[i]-1].x << "/" << verts[faces[i+1]-1].x << "/" << verts[faces[i+2]-1].x << "\n";
		//std::cout << faces[i] << "/" << faces[i+1] << "/" << faces[i+2] << "\n";
		
	}
	//std::cout << "end--------------------------------------------\n";

	//main loop


  glutSwapBuffers();                                      // Draw Frame Buffer
  
}

/***************************************************************************/
void mouse(int button, int state, int x, int y)
{
/* This function I finessed a bit, the value of the printed x,y should
   match the screen, also it remembers where the old value was to avoid multiple
   readings from the same mouse click.  This can cause problems when trying to
   start a line or curve where the last one ended */
        static int oldx = 0;
        static int oldy = 0;
	int mag;
	Color oldColor = selectedColor;

	y *= -1;  //align y with mouse
	y += 500; //ignore 
	mag = (oldx - x)*(oldx - x) + (oldy - y)*(oldy - y);
	if (mag > 20) {
		printf(" x,y is (%d,%d)\n", x,y);

	//if the clicked position is in the wheel, update the color and don't count as a point on the grid
	if (ColorPicker::isToggled)
	{
		selectedColor = cp.updateColor(Point(x,y));

	std::cout << "Color: " << selectedColor.r << " " << selectedColor.g << "  " << selectedColor.b << "\n";
	}	
	if (selectedColor == oldColor)
	{
		selectedColor.randomize();
		//add to point buffer if empty
	//if (listOfPoints.empty())
	//	listOfPoints.emplace_back(Point(x,y));

	//add to point buffer if new point
	
	
	//else
	//	listOfPoints.emplace_back(Point(x,y,selectedColor));
	}
	
	}
	oldx = x;
	oldy = y;
	x_last = x;
	y_last = y;

	

	//std::cout << listOfPoints.size() << " " << listOfPoints[listOfPoints.size()-1].x << " " << listOfPoints[listOfPoints.size()-1].y << "\n";

	
}
 
/***************************************************************************/
void keyboard ( unsigned char key, int x, int y )  // Create Keyboard Function
{

	switch ( key ) {
		case 27:              // When Escape Is Pressed...
			exit ( 0 );   // Exit The Program
			break;        
	        case 'e':             // stub for new screen
			{
				scale = false;
				translation = false;
				rotation = false;
				Transformations::reset();
				glClear(GL_COLOR_BUFFER_BIT);
				verts = vertsOriginal;
				changed = true;
			} break;
			case 'v':
			{
				perspective = !perspective;
				orthographic = !orthographic;
				changed = true; //rerender
				verts = vertsOriginal; //reset points
			} break;
			case 't':
			{
				translation = !translation;
			} break;
			case 'r':
			{
				rotation = !rotation;
			} break;

			case 's':
			{
				if (!translation && !rotation)
				{
					scale = !scale;
				}

				if (rotation) 
				{
					rotateS = true;
					changed = true;
				}
				if (translation) 
				{
					translateS = true;
					changed = true;
				}
				
			} break;

			case 'w':
			{
				if (rotation)
				{ 
					rotateW = true;
					changed = true;
				}
				if (translation) 
				{
					translateW = true;
					changed = true;
				}
			} break;
			case 'a':
			{
				if (scale)
				{
					scaleDown = true;
					changed = true;
				}
				if (rotation)
				{ 
					rotateA = true;
					changed = true;
				}
				if (translation) 
				{
					translateA = true;
					changed = true;
				}
				
			} break;
			case 'd':
			{
				if (scale)
				{
					scaleUp = true;
					changed = true;
				}
				if (rotation)
				{ 
					rotateD = true;
					changed = true;
				}
				if (translation) 
				{
					translateD = true;
					changed = true;
				}
			} break;

		default:       
			break;
	}
}
/***************************************************************************/

int main (int argc, char *argv[])
{
/* This main function sets up the main loop of the program and continues the
   loop until the end of the data is reached.  Then the window can be closed
   using the escape key.						  */

	//read in obj
	if (argc != 2) return 1;
	printf("beginning");
	vertSize = ObjParser::openFile(argv[1]).getVerts().size();
	printf("%d", vertSize);
	
	if (ObjParser::openFile(argv[1]).getFaces().empty())
		std::cout << "empty\n";
	else
	{
		//xfer points to vectors
		for (int i = 0; i < vertSize; ++i)
		{
			verts.emplace_back(vec3<double>(ObjParser::openFile(argv[1]).getVerts()[i].x,
							  				ObjParser::openFile(argv[1]).getVerts()[i].y,
											ObjParser::openFile(argv[1]).getVerts()[i].z
							  ));
			//keep a backup
			vertsOriginal.emplace_back(vec3<double>(ObjParser::openFile(argv[1]).getVerts()[i].x,
							  				ObjParser::openFile(argv[1]).getVerts()[i].y,
											ObjParser::openFile(argv[1]).getVerts()[i].z
							  ));
		}

		for (int i = 0; i < ObjParser::openFile(argv[1]).getFaces().size(); ++i)
		{
			faces.emplace_back(ObjParser::openFile(argv[1]).getFaces()[i].v);
		}
	}       		
	
	glutInit            ( &argc, argv ); 
       	glutInitDisplayMode ( GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH ); 
	glutInitWindowSize  ( WIDTH, HEIGHT ); 
	glutCreateWindow    ( "Computer Graphics" ); 
	glutDisplayFunc     ( display );  
	glutIdleFunc	    ( display );
	glutMouseFunc       ( mouse );
	glutKeyboardFunc    ( keyboard );
        					      
        init_window();				             //create_window

	glutMainLoop        ( );                 // Initialize The Main Loop
}


