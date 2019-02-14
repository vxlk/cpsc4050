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

#define WIDTH 500
#define HEIGHT 500

int x_last,y_last;

//forward declare
void write_pixel(int x, int y, double intensity);

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
	//uint8_t a;
};

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
		if (this->x != rhs.x)
			return true;
		if (this->y != rhs.y)
			return true;
		return false;
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

/*
Connects two points in space
*/
struct Line
{
	static bool isToggled;

	static void draw(const Point& a, Point& b)
	{
		//points to be drawn
		float x = a.x;
		float y = b.y;
		
		//connect old x,y -> new x,y
		//y = (Y2-Y1)/(X2-X1)x + b


		//calculate the number of steps....
		int steps = abs(b.x - a.x) > abs(b.y - a.y) ? 
					abs(b.x - a.x) : 
					abs(b.y - a.y);

		
		//calculate increment values
		float changeX = (b.x - a.x)/(float)steps;
		float changeY = (b.y - a.y)/(float)steps;

		for (int i = 0; i <= steps; ++i)
		{
			write_pixel(x,y,1); //temp intensity COLOR GOES HERE

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

		for (int i = 0; i <= steps; ++i)
		{
			write_pixel(x,y,1); //temp intensity COLOR GOES HERE

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
			//https://buildingvts.com/mathematical-intuition-behind-bezier-curves-2ea4e9645681
			newPoint.x = pow(1.0-t,3.0)*p0.x + 3.0*pow(1.0-t,2.0)*t*p1.x + 3.0*(1-t)*pow(t,2.0)*p2.x + pow(t,3.0)*p3.x;
			newPoint.y = pow(1.0-t,3.0)*p0.y + 3.0*pow(1.0-t,2.0)*t*p1.y + 3.0*(1-t)*pow(t,2.0)*p2.y + pow(t,3.0)*p3.y;
			Line::draw(newPoint, lastPoint);
			lastPoint = newPoint;
		}

	}
};

class ColorPicker
{
	std::vector<Color> palette;

public:
	static bool isToggled;
	ColorPicker() {}
	void draw()
	{
		//draw a circle, each pixel should be a different color
	}

	Color returnColor(Point clickedPos)
	{
		//check the point in the circle, return the color
	}
};

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

void write_pixel(int x, int y, double intensity)
                                         /* Turn on the pixel found at x,y */
{

        //glColor3f (intensity, intensity, intensity);    
		glColor3f(.5, .3, .1);             
        glBegin(GL_POINTS);
           glVertex3i( x, y, 0);
        glEnd();	
}

//***************************************************************************/

void display ( void )   // Create The Display Function
{

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	      // Clear Screen 

  //draw a pixel at mouse pos
  write_pixel(x_last,y_last,1.0);

  // CALL YOUR CODE HERE
  int pos = listOfPoints.size()-1;
  static int lastPos = 4;
  if (listOfPoints.size() >= 2 && Line::isToggled)
		for (; pos > 0; --pos)
			Line::draw(pos);
  if (listOfPoints.size() >= 4  && BezierCurve::isToggled)
  {
	
	if (listOfPoints.size()%4==0)
		if (pos != lastPos)
			lastPos = listOfPoints.size();
	
	for(int i = lastPos; i >= 4; i-=4)
		BezierCurve::draw(i-1);
  }
  
  /*
  {
		//for (; pos > 3; --pos)
		if (listOfPoints.size()%4==0)
			if (indexes.empty() || pos != indexes[indexes.size()-1])
				indexes.push_back(pos);
		for (unsigned long int i = 0; i < indexes.size(); ++i)
			BezierCurve::draw(indexes[i]);

		//std::cout << indexes.size() << "\n";
  }
  */

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
	
	y *= -1;  //align y with mouse
	y += 500; //ignore 
	mag = (oldx - x)*(oldx - x) + (oldy - y)*(oldy - y);
	if (mag > 20) {
		printf(" x,y is (%d,%d)\n", x,y);
	}
	oldx = x;
	oldy = y;
	x_last = x;
	y_last = y;

	//add to point buffer if empty
	if (listOfPoints.empty())
		listOfPoints.emplace_back(Point(x,y));

	//add to point buffer if new point
	if (Point(x,y) != listOfPoints[listOfPoints.size()-1])
		listOfPoints.emplace_back(Point(x,y));
	
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
		        glClearColor(0.0,0.0,0.0,0.0);
				glClear(GL_COLOR_BUFFER_BIT);
			} break;
			case 'l':
			{
				Line::isToggled = true;
				if (BezierCurve::isToggled)
					BezierCurve::isToggled = false;
			} break;
			case 'c':
			{
				BezierCurve::isToggled = true;
				if (Line::isToggled)
					Line::isToggled = false;
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
	
	glutInit            ( &argc, argv ); 
       	glutInitDisplayMode ( GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH ); 
	glutInitWindowSize  ( 500,500 ); 
	glutCreateWindow    ( "Computer Graphics" ); 
	glutDisplayFunc     ( display );  
	glutIdleFunc	    ( display );
	glutMouseFunc       ( mouse );
	glutKeyboardFunc    ( keyboard );
        					      
        init_window();				             //create_window
						       		
	glutMainLoop        ( );                 // Initialize The Main Loop
}


