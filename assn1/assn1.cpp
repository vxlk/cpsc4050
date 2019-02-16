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

#define WIDTH 500
#define HEIGHT 500

int x_last,y_last;

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

	static void draw(const Point& a, const Point& b)
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

		Color incrementer = b.c - a.c;
		incrementer = incrementer /(double) steps;
		//calculate increment values
		float changeX = (b.x - a.x)/(float)steps;
		float changeY = (b.y - a.y)/(float)steps;

		for (int i = 0; i <= steps; ++i)
		{
			write_pixel(x,y,b.c+incrementer); //temp intensity COLOR GOES HERE

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
			write_pixel(x,y,selectedColor); //temp intensity COLOR GOES HERE

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
		bool useRed = false;
		bool useGreen = false;
		bool useBlue = false;
		double r = 0.0;
		double g = 0.0;
		double b = 0.0;

		/*
		//draw a circle, each pixel should be a different color
		for(int i = 0; i < 255; ++i)
			for(int u = 0; u < 255; ++u)
			{
			if (r < 1.0 && !useRed) r+=inc;
			else if (g < 1.0 && !useGreen)
			{ 
				g+=inc;
				r = 0;
				useRed = true;
			}
			else if (b < 1.0 && !useBlue)
			{
				b += inc;
				g = 0;
				useGreen = true;
			}

			else if (useRed && useGreen && r < 1.0)
				while(r != 1.0)
				{
					g += inc;

				}

			//u wanna make this a pair operator < > point
			*/
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
  write_pixel(x_last,y_last,selectedColor);

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
  
  if (ColorPicker::isToggled)
  	cp.draw();

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
	if (listOfPoints.empty())
		listOfPoints.emplace_back(Point(x,y));

	//add to point buffer if new point
	
	//if (Point(x,y) != listOfPoints[listOfPoints.size()-1])
	else
		listOfPoints.emplace_back(Point(x,y,selectedColor));
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
				glClear(GL_COLOR_BUFFER_BIT);
				listOfPoints.clear();
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
			case 'r':
				ColorPicker::isToggled = !ColorPicker::isToggled;
			break;
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


