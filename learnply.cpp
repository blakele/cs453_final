#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <fstream>
#include <vector>

#define _USE_MATH_DEFINES

#include "glError.h"
#include "gl/glew.h"
#include "gl/freeglut.h"
#include "ply.h"
#include "icVector.H"
#include "icMatrix.H"
#include "polyhedron.h"
#include "polyline.h"
#include "trackball.h"
#include "tmatrix.h"
#include <cmath>
#include <iostream>
#include <tuple>

using namespace std;

Polyhedron* poly;

bool showPickedPoint;
icVector3 pickedPoint;

//transfer matrix used to get xyz
double matrix[3][3] = { {0.4124564, 0.2126729, 0.0193339},
								{0.3575761, 0.7151522, 0.1191920},
								{0.1804375, 0.0721750, 0.9503041} };

double inverse_matrix[3][3] = { { 3.24045484, -0.96926639,  0.05564342},
								{-1.53713885,  1.87601093, -0.20402585},
								{-0.49853155,  0.04155608,  1.05722516} };
//used for converting to xyz2lab


/*scene related variables*/
const float zoomspeed = 0.9;
const int view_mode = 0;		// 0 = othogonal, 1=perspective
const double radius_factor = 1.0;
int win_width = 800;
int win_height = 800;
float aspectRatio = win_width / win_height;
/*
Use keys 1 to 0 to switch among different display modes.
Each display mode can be designed to show one type 
visualization result.

Predefined ones: 
display mode 1: solid rendering
display mode 2: show wireframes
display mode 3: render each quad with colors of vertices
*/
int display_mode = 1;

/*User Interaction related variabes*/
float s_old, t_old;
float rotmat[4][4];
double zoom = 1.0;
double translation[2] = { 0, 0 };
int mouse_mode = -2;	// -1 = no action, 1 = tranlate y, 2 = rotate

/*IBFV related variables*/
//https://www.win.tue.nl/~vanwijk/ibfv/
#define	NPN 64
#define SCALE 4.0
int    Npat = 32;
int    iframe = 0;
float  tmax = win_width / (SCALE*NPN);
float  dmax = SCALE / win_width;
unsigned char *pixels;

#define DM  ((float) (1.0/(100-1.0)))

/* Contour Variables */

std::vector<std::pair<PolyLine, float>> contours;
const int num_intervals = 20;

/******************************************************************************
Forward declaration of functions
******************************************************************************/

void init(void);
void makePatterns(void);

/*glut attaching functions*/
void keyboard(unsigned char key, int x, int y);
void motion(int x, int y);
void display(void);
void mouse(int button, int state, int x, int y);
void mousewheel(int wheel, int direction, int x, int y);
void reshape(int width, int height);

/*functions for element picking*/
void display_vertices(GLenum mode, Polyhedron* poly);
void display_quads(GLenum mode, Polyhedron* poly);
void display_selected_vertex(Polyhedron* poly);
void display_selected_quad(Polyhedron* poly);

/*display vis results*/
void display_polyhedron(Polyhedron* poly);

//get linear rgb values from sRGB
double* linear_rgb(double r, double g, double b) {
	double rgb[3] = { r, g, b };
	double* linearRgb = new double[3];
	for (int i = 0; i < 3; i++) {
		double value = rgb[i] / 255.;
		if (value > 0.04045) {
			value = pow(((value + 0.055) / 1.055), 2.4);
		}
		else {
			value = value / 12.92;
		}
		linearRgb[i] = value * 100.;
	}
	return linearRgb;
}

double* linear_rgb(Vertex* v) {
	return linear_rgb(v->R, v->G, v->B);

}

double* sRGB(double r, double g, double b) {
	double linearRgb[3] = { r, g, b };
	double* sRgb = new double[3];
	for (int i = 0; i < 3; i++) {
		double value = linearRgb[i] / 100.;
		if (value > 0.00313080495356037152) {
			value = (1.055 * pow(value, 1. / 2.4)) - 0.055;
		}
		else {
			value = value * 12.92;
		}
		sRgb[i] = value * 255.;
	}
	return sRgb;
}

//this function returns the vertex specified by user input
Vertex* get_vertex(double x, double y) {
	Quad* q = NULL;
	double x1, y1, x2, y2;
	for (int i = 0; i < poly->nquads; i++) {
		q = poly->qlist[i];
		x1 = q->verts[0]->x;
		y1 = q->verts[0]->y;
		x2 = q->verts[0]->x;
		y2 = q->verts[0]->y;

		for (int j = 0; j < 4; j++) {
			if (x1 >= q->verts[j]->x) {
				x1 = q->verts[j]->x;
			}
			if (x2 <= q->verts[j]->x) {
				x2 = q->verts[j]->x;
			}
			if (y1 >= q->verts[j]->y) {
				y1 = q->verts[j]->y;
			}
			if (y2 <= q->verts[j]->y) {
				y2 = q->verts[j]->y;
			}
		}
		//similar to how we find x1,x2,y1,y2 in hw2
		if (x <= x2 && x >= x1 && y <= y2 && y >= y1) {
			break;
		}
	}

	for (int i = 0; i < 4; i++) {
		//cout << "vertex " << i << endl;
		//cout << "x = " << q->verts[i]->x << ", y = " << q->verts[i]->y << endl;
		if (x1 == q->verts[i]->x && y1 == q->verts[i]->y) {
			return q->verts[i];
		}
	}
}

//this function prints the rgb value at a specific vertex
void print_rgb(Vertex* v) {
	printf("R: %f, G: %f, B: %f\n", v->R, v->G, v->B);
}


double* rgb2xyz(double r, double g, double b) {
	double temp;
	double* xyz_array = new double[3];
	for (int i = 0; i < 3; i++) {
		temp = ((r * matrix[0][i]) + (g * matrix[1][i]) + (b * matrix[2][i]));
		xyz_array[i] = temp;
	}

	// column/row backwards for dot product
	/*for (int i = 0; i < 3; i++) {
		temp = ((r * matrix[i][0]) + (g * matrix[i][1]) + (b * matrix[i][2]));
		xyz_array[i] = temp;
	}*/

	//printf("printing xyz\n");
	//printf("X: %f, Y: %f, Z: %f\n", xyz_array[0], xyz_array[1], xyz_array[2]);

	return xyz_array;
}

//get xyz space from linear rgb values
double* rgb2xyz(Vertex* v) {
	//initialize variables
	double a, b, c;
	a = v->R;
	b = v->G;
	c = v->B;
	
	return rgb2xyz(a, b, c);
}

//helper function for xyz2lab function
double f(double x) {
	double limit = 0.008856;
	if (x > limit) {
		return pow(x, (1. / 3.));
	}
	else {
		return 7.787 * x + 16. / 116.;
	}
}

double finverse(double x) {
	double limit = 0.008856;
	double a = 7.787;
	double b = 16. / 116.;
	double ylim = a * limit + b;
	if (x > ylim) {
		return x * x * x;
	}
	else {
		return (x - b) / a;
	}
}


double* xyz2lab(double x, double y, double z) {
	double L, A, B;
	double xn = 95.047;
	double yn = 100.0;
	double zn = 108.883;

	L = 116. * (f(y / yn) - (16. / 116.));
	A = 500. * (f(x / xn) - f(y / yn));
	B = 200. * (f(y / yn) - f(z / zn));

	double* lab = new double[3];
	lab[0] = L;
	lab[1] = A;
	lab[2] = B;

	//printf("printing LAB values for vertex\n");
	//printf("L: %f, A: %f, B: %f\n", lab[0], lab[1], lab[2]);
	return lab;
}

//this function takes a vertex and returns the LAB values from the vertex's rgb
double* cielab(Vertex* v) {
	double* temp = rgb2xyz(v);
	return xyz2lab(temp[0], temp[1], temp[2]);
	delete [] temp;
}

double* lab2msh(double L, double A, double B) {
	double temp = (L * L) + (A * A) + (B * B);

	double* msh = new double[3];
	msh[0] = sqrt(temp);
	msh[1] = acos(L / msh[0]);
	msh[2] = atan(B / A);
	//printf("printing msh values for vertex\n");
	//printf("M: %f, S: %f, H: %f\n", msh[0], msh[1], msh[2]);
	return msh;
}

double* rgb2msh(double r, double g, double b) {
	double* xyz = rgb2xyz(r, g, b);
	double* lab = xyz2lab(xyz[0], xyz[1], xyz[2]);
	delete[] xyz;
	double L = lab[0];
	double A = lab[1];
	double B = lab[2];
	delete[] lab;
	double* msh = lab2msh(L, A, B);
	return msh;
}

//this function uses LAB variables to get rgb2msh
double* rgb2msh(Vertex* v) {
	return rgb2msh(v->R, v->G, v->B);
}

double* xyz2rgb(double x, double y, double z) {
	double temp;
	double* rgb_array = new double[3];
	for (int i = 0; i < 3; i++) {
		temp = ((x * inverse_matrix[0][i]) + (y * inverse_matrix[1][i]) + (z * inverse_matrix[2][i]));
		rgb_array[i] = temp;
	}
	return rgb_array;
}

double* lab2xyz(double l, double a, double b) {
	double xn = 95.047;
	double yn = 100.0;
	double zn = 108.883;

	double* xyz_array = new double[3];
	xyz_array[0] = xn * finverse((a / 500.) + (l + 16.) / 116.);
	xyz_array[1] = yn * finverse((l + 16.) / 116.);
	xyz_array[2] = zn * finverse((l + 16.) / 116. - (b / 200.));

	return xyz_array;
}


double* msh2lab(double m, double s, double h) {
	double* lab = new double[3];
	lab[0] = m * cos(s);
	lab[1] = m * sin(s) * cos(h);
	lab[2] = m * sin(s) * sin(h);
	return lab;
}

double* msh2rgb(double m, double s, double h) {
	double* lab_arr = msh2lab(m, s, h);
	double* xyz_arr = lab2xyz(lab_arr[0], lab_arr[1], lab_arr[2]);
	delete[] lab_arr;
	double* rgb_arr = xyz2rgb(xyz_arr[0], xyz_arr[1], xyz_arr[2]);
	delete[] xyz_arr;
	return rgb_arr;
}

double hueSpin(double Msat, double ssat, double hsat, double Munsat) {
	if (Msat >= Munsat) {
		return hsat;
	}
	else {
		double hSpin = ssat * sqrt((Munsat * Munsat) - (Msat * Msat)) / (Msat * sin(ssat));
		if (hsat > -M_PI/3.) {
			return hsat + hSpin;
		}
		else {
			return hsat - hSpin;
		}
	}
}

double* interpolateColor(double* RGB1, double* RGB2, double interp) {
	double* linearRGB1 = linear_rgb(RGB1[0], RGB1[1], RGB1[2]);
	double* msh1 = rgb2msh(linearRGB1[0], linearRGB1[1], linearRGB1[2]);
	double M1 = msh1[0];
	double s1 = msh1[1];
	double h1 = msh1[2];
	delete[] linearRGB1;
	delete[] msh1;

	double* linearRGB2 = linear_rgb(RGB2[0], RGB2[1], RGB2[2]);
	double* msh2 = rgb2msh(linearRGB2[0], linearRGB2[1], linearRGB2[2]);
	double M2 = msh2[0];
	double s2 = msh2[1];
	double h2 = msh2[2];
	delete[] linearRGB2;
	delete[] msh2;

	// two saturated points distinct in color, add white between them
	//printf("h diff: %f", fabs(h1 - h2));
	if (s1 > 0.05 && s2 > 0.05 && fabs(h1 - h2) > (M_PI / 3.)) {
		double Mmid = M1;
		if (M2 > Mmid) Mmid = M2;
		if (88. > Mmid) Mmid = 88.;
		if (interp < 0.5) {
			M2 = Mmid;
			s2 = 0.;
			h2 = 0;
			interp = 2 * interp;
		} 
		else {
			M1 = Mmid;
			s1 = 0.;
			h1 = 0;
			interp = 2 * interp - 1;
		}
	}
	if (s1 < 0.05 && s2 > 0.05) {
		h1 = hueSpin(M2, s2, h2, M1);
	}
	else if (s2 < 0.05 && s1 > 0.05) {
		h2 = hueSpin(M1, s1, h1, M2);
	}

	double MshMid[3];
	MshMid[0] = (1 - interp) * M1 + interp * M2;
	MshMid[1] = (1 - interp) * s1 + interp * s2;
	MshMid[2] = (1 - interp) * h1 + interp * h2;

	double* rgb_arr = msh2rgb(MshMid[0], MshMid[1], MshMid[2]);
	double* srgb_arr = sRGB(rgb_arr[0], rgb_arr[1], rgb_arr[2]);
	delete[] rgb_arr;
	return srgb_arr;
}

void get_contours(int height_mul) {
	contours.clear();
	// find out max and min scalar values
	double min = poly->vlist[0]->scalar;
	double max = poly->vlist[0]->scalar;
	for (int i = 0; i < poly->nverts; i++) {
		Vertex* temp_v = poly->vlist[i];
		double scalar = temp_v->scalar;
		if (scalar > max) {
			max = scalar;
		}
		else if (scalar < min) {
			min = scalar;
		}
	}

	// iterate through once for each interval, and create a contour for that s_0
	for (int j = 0; j < num_intervals; j++) {
		double s_0 = (((double)j / (double)num_intervals) * (max - min)) + min;
		PolyLine contour;
		printf("s_0 = %f\n", s_0);

		// iterate through each vertex and classify them by if the are bigger or smaller than s_0 
		for (int i = 0; i < poly->nverts; i++) {
			Vertex* v = poly->vlist[i];
			if (v->scalar > s_0) {
				v->vertex_type = 1;
			}
			else {
				v->vertex_type = -1;
			}
			v->z = (v->scalar - min) / (max - min) * height_mul;
		}

		// iterate through each edge, and find the ones with different vertex types on each side
		// then interpolate to find the intersection and create a vertex there
		for (int i = 0; i < poly->nedges; i++) {
			Edge* temp_e = poly->elist[i];
			if (temp_e->verts[0]->vertex_type != temp_e->verts[1]->vertex_type) {
				double alpha = (s_0 - temp_e->verts[0]->scalar) / (temp_e->verts[1]->scalar - temp_e->verts[0]->scalar);
				double v_x = alpha * (temp_e->verts[1]->x - temp_e->verts[0]->x) + temp_e->verts[0]->x;
				double v_y = alpha * (temp_e->verts[1]->y - temp_e->verts[0]->y) + temp_e->verts[0]->y;
				double v_z = alpha * (temp_e->verts[1]->z - temp_e->verts[0]->z) + temp_e->verts[0]->z;
				Vertex* v = new Vertex(v_x, v_y, v_z);
				temp_e->is_crossing = 1;
				temp_e->crossing = v;
			}
		}

		// iterate through all the quads and their edges, creating lines between any edges in a quad with crossing vertexes
		for (int i = 0; i < poly->nquads; i++) {
			Quad* q = poly->qlist[i];
			Vertex* last_vertex = NULL;

			for (int k = 0; k < 4; k++) {
				Edge* e = q->edges[k];
				if (e->is_crossing > 0) {
					if (last_vertex) {
						LineSegment line(e->crossing->x, e->crossing->y, e->crossing->z, last_vertex->x, last_vertex->y, last_vertex->z);
						(contour).push_back(line);
					}
					last_vertex = e->crossing;
				}
			}
		}

		// remove the vertex_type and crossing vertexes for this s_0 contour so we can find the next countour
		for (int i = 0; i < poly->nquads; i++) {
			Quad* q = poly->qlist[i];

			for (int k = 0; k < 4; k++) {
				Edge* e = q->edges[k];
				if (e->crossing)
					delete(e->crossing);
				e->crossing = NULL;
				e->is_crossing = 0;
				e->verts[0]->vertex_type = 0;
				e->verts[1]->vertex_type = 0;
			}
		}

		contours.push_back(std::pair<PolyLine, float>(contour, s_0));
	}
}

/*display utilities*/

/*
draw a sphere
x, y, z are the coordiate of the dot
radius of the sphere 
R: the red channel of the color, ranges [0, 1]
G: the green channel of the color, ranges [0, 1]
B: the blue channel of the color, ranges [0, 1]
*/
void drawDot(double x, double y, double z, double radius = 0.1, double R = 1.0, double G = 0.0, double B = 0.0) {

	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1., 1.);
	glEnable(GL_DEPTH_TEST);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glShadeModel(GL_SMOOTH);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);

	GLfloat mat_diffuse[4];

	{
		mat_diffuse[0] = R;
		mat_diffuse[1] = G;
		mat_diffuse[2] = B;
		mat_diffuse[3] = 1.0;
	}

	glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);

	GLUquadric* quad = gluNewQuadric();

	glPushMatrix();
	glTranslatef(x, y, z);
	gluSphere(quad, radius, 50, 50);
	glPopMatrix();

	gluDeleteQuadric(quad);
}

/*
draw a line segment
width: the width of the line, should bigger than 0
R: the red channel of the color, ranges [0, 1]
G: the green channel of the color, ranges [0, 1]
B: the blue channel of the color, ranges [0, 1]
*/
void drawLineSegment(LineSegment ls, double width = 1.0, double R = 1.0, double G = 0.0, double B = 0.0) {

	glDisable(GL_LIGHTING);
	glEnable(GL_LINE_SMOOTH);
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glLineWidth(width);

	glBegin(GL_LINES);
	glColor3f(R, G, B);
	glVertex3f(ls.start.x, ls.start.y, ls.start.z);
	glVertex3f(ls.end.x, ls.end.y, ls.end.z);
	glEnd();

	glDisable(GL_BLEND);
}

/*
draw a polyline
width: the width of the line, should bigger than 0
R: the red channel of the color, ranges [0, 1]
G: the green channel of the color, ranges [0, 1]
B: the blue channel of the color, ranges [0, 1]
*/
void drawPolyline(PolyLine pl, double width = 1.0, double R = 1.0, double G = 0.0, double B = 0.0) {
	
	glDisable(GL_LIGHTING);
	glEnable(GL_LINE_SMOOTH);
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glLineWidth(width);

	glBegin(GL_LINES);
	glColor3f(R, G, B);

	for (int i = 0; i < pl.size(); i++) {
		glVertex3f(pl[i].start.x, pl[i].start.y, pl[i].start.z);
		glVertex3f(pl[i].end.x, pl[i].end.y, pl[i].end.z);
	}

	glEnd();

	glDisable(GL_BLEND);
}

/******************************************************************************
Main program.
******************************************************************************/
int main(int argc, char* argv[])
{
	/*load mesh from ply file*/
	//FILE* this_file = fopen("../quadmesh_2D/vector_data/saddle.ply", "r");
	FILE* this_file = fopen("../quadmesh_2D/scalar_data/2x_square_plus_y_square.ply", "r");

	poly = new Polyhedron(this_file);
	fclose(this_file);
	
	/*initialize the mesh*/
	poly->initialize(); // initialize the mesh
	poly->write_info();


	/*init glut and create window*/
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
	glutInitWindowPosition(20, 20);
	glutInitWindowSize(win_width, win_height);
	glutCreateWindow("Scientific Visualization");


	/*initialize openGL*/
	init();

	/*prepare the noise texture for IBFV*/
	makePatterns();
	
	/*the render function and callback registration*/
	glutKeyboardFunc(keyboard);
	glutReshapeFunc(reshape);
	glutDisplayFunc(display);
	glutIdleFunc(display);
	glutMotionFunc(motion);
	glutMouseFunc(mouse);
	glutMouseWheelFunc(mousewheel);
	
	/*event processing loop*/
	glutMainLoop();
	
	/*clear memory before exit*/
	poly->finalize();	// finalize everything
	free(pixels);
	return 0;
}


/******************************************************************************
Set projection mode
******************************************************************************/

void set_view(GLenum mode)
{
	GLfloat light_ambient0[] = { 0.3, 0.3, 0.3, 1.0 };
	GLfloat light_diffuse0[] = { 0.7, 0.7, 0.7, 1.0 };
	GLfloat light_specular0[] = { 0.0, 0.0, 0.0, 1.0 };

	GLfloat light_ambient1[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat light_diffuse1[] = { 0.5, 0.5, 0.5, 1.0 };
	GLfloat light_specular1[] = { 0.0, 0.0, 0.0, 1.0 };

	glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient0);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse0);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular0);

	glLightfv(GL_LIGHT1, GL_AMBIENT, light_ambient1);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, light_diffuse1);
	glLightfv(GL_LIGHT1, GL_SPECULAR, light_specular1);


	glMatrixMode(GL_PROJECTION);
	if (mode == GL_RENDER)
		glLoadIdentity();

	if (aspectRatio >= 1.0) {
		if (view_mode == 0)
			glOrtho(-radius_factor * zoom * aspectRatio, radius_factor * zoom * aspectRatio, -radius_factor * zoom, radius_factor * zoom, -1000, 1000);
		else
			glFrustum(-radius_factor * zoom * aspectRatio, radius_factor * zoom * aspectRatio, -radius_factor* zoom, radius_factor* zoom, 0.1, 1000);
	}
	else {
		if (view_mode == 0)
			glOrtho(-radius_factor * zoom, radius_factor * zoom, -radius_factor * zoom / aspectRatio, radius_factor * zoom / aspectRatio, -1000, 1000);
		else
			glFrustum(-radius_factor * zoom, radius_factor * zoom, -radius_factor* zoom / aspectRatio, radius_factor* zoom / aspectRatio, 0.1, 1000);
	}


	GLfloat light_position[3];
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	light_position[0] = 5.5;
	light_position[1] = 0.0;
	light_position[2] = 0.0;
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
	light_position[0] = -0.1;
	light_position[1] = 0.0;
	light_position[2] = 0.0;
	glLightfv(GL_LIGHT2, GL_POSITION, light_position);
}

/******************************************************************************
Update the scene
******************************************************************************/

void set_scene(GLenum mode, Polyhedron* poly)
{
	glTranslatef(translation[0], translation[1], -3.0);

	/*multiply rotmat to current mat*/
	{
		int i, j, index = 0;

		GLfloat mat[16];

		for (i = 0; i < 4; i++)
			for (j = 0; j < 4; j++)
				mat[index++] = rotmat[i][j];

		glMultMatrixf(mat);
	}

	glScalef(0.9 / poly->radius, 0.9 / poly->radius, 0.9 / poly->radius);
	glTranslatef(-poly->center.entry[0], -poly->center.entry[1], -poly->center.entry[2]);
}


/******************************************************************************
Init scene
******************************************************************************/

void init(void) {

	mat_ident(rotmat);

	/* select clearing color */
	glClearColor(0.0, 0.0, 0.0, 0.0);  // background
	glShadeModel(GL_FLAT);
	glPolygonMode(GL_FRONT, GL_FILL);

	glDisable(GL_DITHER);
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);
	
	//set pixel storage modes
	glPixelStorei(GL_PACK_ALIGNMENT, 1);
	
	glEnable(GL_NORMALIZE);
	if (poly->orientation == 0)
		glFrontFace(GL_CW);
	else
		glFrontFace(GL_CCW);
}


/******************************************************************************
Pick objects from the scene
******************************************************************************/

int processHits(GLint hits, GLuint buffer[])
{
	unsigned int i, j;
	GLuint names, * ptr;
	double smallest_depth = 1.0e+20, current_depth;
	int seed_id = -1;
	unsigned char need_to_update;

	ptr = (GLuint*)buffer;
	for (i = 0; i < hits; i++) {  /* for each hit  */
		need_to_update = 0;
		names = *ptr;
		ptr++;

		current_depth = (double)*ptr / 0x7fffffff;
		if (current_depth < smallest_depth) {
			smallest_depth = current_depth;
			need_to_update = 1;
		}
		ptr++;
		current_depth = (double)*ptr / 0x7fffffff;
		if (current_depth < smallest_depth) {
			smallest_depth = current_depth;
			need_to_update = 1;
		}
		ptr++;
		for (j = 0; j < names; j++) {  /* for each name */
			if (need_to_update == 1)
				seed_id = *ptr - 1;
			ptr++;
		}
	}
	return seed_id;
}

/******************************************************************************
Diaplay all quads for selection
******************************************************************************/

void display_quads(GLenum mode, Polyhedron* this_poly)
{
	unsigned int i, j;
	GLfloat mat_diffuse[4];

	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1., 1.);
	glEnable(GL_DEPTH_TEST);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glShadeModel(GL_SMOOTH);
	//glDisable(GL_LIGHTING);

	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

	for (i = 0; i < this_poly->nquads; i++) {
		if (mode == GL_SELECT)
			glLoadName(i + 1);

		Quad* temp_q = this_poly->qlist[i];
		{
			mat_diffuse[0] = 1.0;
			mat_diffuse[1] = 1.0;
			mat_diffuse[2] = 0.0;
			mat_diffuse[3] = 1.0;
		}
		glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
		
		glBegin(GL_POLYGON);
		for (j = 0; j < 4; j++) {
			Vertex* temp_v = temp_q->verts[j];
			//glColor3f(0, 0, 0);
			glVertex3d(temp_v->x, temp_v->y, temp_v->z);
		}
		glEnd();
	}
}

/******************************************************************************
Diaplay all vertices for selection
******************************************************************************/

void display_vertices(GLenum mode, Polyhedron* this_poly)
{
	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1., 1.);
	glEnable(GL_DEPTH_TEST);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glShadeModel(GL_SMOOTH);
	glDisable(GL_LIGHTING);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

	for (int i = 0; i < this_poly->nverts; i++) {
		if (mode == GL_SELECT)
			glLoadName(i + 1);

		Vertex* temp_v = this_poly->vlist[i];

		{
			GLUquadric* quad = gluNewQuadric();

			glPushMatrix();
			glTranslatef(temp_v->x, temp_v->y, temp_v->z);
			glColor4f(0, 0, 1, 1.0);
			gluSphere(quad, this_poly->radius * 0.01, 50, 50);
			glPopMatrix();

			gluDeleteQuadric(quad);
		}
	}
}

/******************************************************************************
Diaplay selected quad
******************************************************************************/

void display_selected_quad(Polyhedron* this_poly)
{
	if (this_poly->selected_quad == -1)
	{
		return;
	}

	unsigned int i, j;
	GLfloat mat_diffuse[4];

	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1., 1.);
	glDisable(GL_DEPTH_TEST);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glShadeModel(GL_SMOOTH);
	glDisable(GL_LIGHTING);

	Quad* temp_q = this_poly->qlist[this_poly->selected_quad];

	glBegin(GL_POLYGON);
	for (j = 0; j < 4; j++) {
		Vertex* temp_v = temp_q->verts[j];
		glColor3f(1.0, 0.0, 1.0);
		glVertex3d(temp_v->x, temp_v->y, 0.0);
	}
	glEnd();
}

/******************************************************************************
Diaplay selected vertex
******************************************************************************/

void display_selected_vertex(Polyhedron* this_poly)
{
	if (this_poly->selected_vertex == -1)
	{
		return;
	}

	Vertex* temp_v = this_poly->vlist[this_poly->selected_vertex];

	drawDot(temp_v->x, temp_v->y, temp_v->z, this_poly->radius * 0.01, 1.0, 0.0, 0.0);

}


/******************************************************************************
Callback function for glut window reshaped
******************************************************************************/

void reshape(int width, int height) {

	win_width = width;
	win_height = height;

	aspectRatio = (float)width / (float)height;

	glViewport(0, 0, width, height);

	set_view(GL_RENDER);

	/*Update pixels buffer*/
	free(pixels);
	pixels = (unsigned char *)malloc(sizeof(unsigned char)*win_width*win_height * 3);
	memset(pixels, 255, sizeof(unsigned char)*win_width*win_height * 3);
}


/******************************************************************************
Callback function for dragging mouse
******************************************************************************/

void motion(int x, int y) {
	float r[4];
	float s, t;

	s = (2.0 * x - win_width) / win_width;
	t = (2.0 * (win_height - y) - win_height) / win_height;

	if ((s == s_old) && (t == t_old))
		return;

	switch (mouse_mode) {
	case 2:

		Quaternion rvec;

		mat_to_quat(rotmat, rvec);
		trackball(r, s_old, t_old, s, t);
		add_quats(r, rvec, rvec);
		quat_to_mat(rvec, rotmat);

		s_old = s;
		t_old = t;

		display();
		break;

	case 1:

		translation[0] += (s - s_old);
		translation[1] += (t - t_old);

		s_old = s;
		t_old = t;

		display();
		break;
	}
}

/******************************************************************************
Callback function for mouse clicks
******************************************************************************/

void mouse(int button, int state, int x, int y) {

	int key = glutGetModifiers();

	if (button == GLUT_LEFT_BUTTON || button == GLUT_RIGHT_BUTTON) {
		
		if (state == GLUT_DOWN) {
			float xsize = (float)win_width;
			float ysize = (float)win_height;

			float s = (2.0 * x - win_width) / win_width;
			float t = (2.0 * (win_height - y) - win_height) / win_height;

			s_old = s;
			t_old = t;

			/*translate*/
			if (button == GLUT_LEFT_BUTTON)
			{
				mouse_mode = 1;
			}

			/*rotate*/
			if (button == GLUT_RIGHT_BUTTON)
			{
				mouse_mode = 2;
			}
		}
		else if (state == GLUT_UP) {

			if (button == GLUT_LEFT_BUTTON && key == GLUT_ACTIVE_SHIFT) {  // build up the selection feedback mode

				/*select face*/

				GLuint selectBuf[512];
				GLint hits;
				GLint viewport[4];

				glGetIntegerv(GL_VIEWPORT, viewport);

				glSelectBuffer(win_width, selectBuf);
				(void)glRenderMode(GL_SELECT);

				glInitNames();
				glPushName(0);

				glMatrixMode(GL_PROJECTION);
				glPushMatrix();
				glLoadIdentity();

				/*create 5x5 pixel picking region near cursor location */
				gluPickMatrix((GLdouble)x, (GLdouble)(viewport[3] - y), 1.0, 1.0, viewport);

				set_view(GL_SELECT);
				set_scene(GL_SELECT, poly);
				display_quads(GL_SELECT, poly);

				glMatrixMode(GL_PROJECTION);
				glPopMatrix();
				glFlush();

				glMatrixMode(GL_MODELVIEW);

				hits = glRenderMode(GL_RENDER);
				poly->selected_quad = processHits(hits, selectBuf);
				printf("Selected quad id = %d\n", poly->selected_quad);
				glutPostRedisplay();

			}
			else if (button == GLUT_LEFT_BUTTON && key == GLUT_ACTIVE_CTRL)
			{
				/*select vertex*/

				GLuint selectBuf[512];
				GLint hits;
				GLint viewport[4];

				glGetIntegerv(GL_VIEWPORT, viewport);

				glSelectBuffer(win_width, selectBuf);
				(void)glRenderMode(GL_SELECT);

				glInitNames();
				glPushName(0);

				glMatrixMode(GL_PROJECTION);
				glPushMatrix();
				glLoadIdentity();

				/*  create 5x5 pixel picking region near cursor location */
				gluPickMatrix((GLdouble)x, (GLdouble)(viewport[3] - y), 1.0, 1.0, viewport);

				set_view(GL_SELECT);
				set_scene(GL_SELECT, poly);
				display_vertices(GL_SELECT, poly);

				glMatrixMode(GL_PROJECTION);
				glPopMatrix();
				glFlush();

				glMatrixMode(GL_MODELVIEW);

				hits = glRenderMode(GL_RENDER);
				poly->selected_vertex = processHits(hits, selectBuf);
				printf("Selected vert id = %d\n", poly->selected_vertex);
				glutPostRedisplay();

			}

			mouse_mode = -1;
		}
	}
}

/******************************************************************************
Callback function for mouse wheel scroll
******************************************************************************/

void mousewheel(int wheel, int direction, int x, int y) {
	if (direction == 1) {
		zoom *= zoomspeed;
		glutPostRedisplay();
	}
	else if (direction == -1) {
		zoom /= zoomspeed;
		glutPostRedisplay();
	}
}

/*Display IBFV*/
void makePatterns(void)
{
	pixels = (unsigned char *)malloc(sizeof(unsigned char)*win_width*win_height * 3);
	memset(pixels, 255, sizeof(unsigned char)*win_width*win_height * 3);

	int lut[256];
	int phase[NPN][NPN];
	GLubyte pat[NPN][NPN][4];
	int i, j, k, t;

	for (i = 0; i < 256; i++) lut[i] = i < 127 ? 0 : 255;
	for (i = 0; i < NPN; i++)
		for (j = 0; j < NPN; j++) phase[i][j] = rand() % 256;

	for (k = 0; k < Npat; k++) {
		t = k * 256 / Npat;
		for (i = 0; i < NPN; i++)
			for (j = 0; j < NPN; j++) {
				pat[i][j][0] =
					pat[i][j][1] =
					pat[i][j][2] = lut[(t + phase[i][j]) % 255];
				pat[i][j][3] = int(0.12 * 255);
			}
		glNewList(k + 1, GL_COMPILE);
		glTexImage2D(GL_TEXTURE_2D, 0, 4, NPN, NPN, 0, GL_RGBA, GL_UNSIGNED_BYTE, pat);
		glEndList();
	}

}

void displayIBFV(void)
{
	glDisable(GL_LIGHTING);
	glDisable(GL_LIGHT0);
	glDisable(GL_LIGHT1);
	glDisable(GL_POLYGON_OFFSET_FILL);
	glDisable(GL_DEPTH_TEST);

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);

	glEnable(GL_TEXTURE_2D);
	glShadeModel(GL_FLAT);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glClearColor(1.0, 1.0, 1.0, 1.0);  // background for rendering color coding and lighting
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	/*draw the model with using the pixels, using vector field to advert the texture coordinates*/
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, win_width, win_height, 0, GL_RGB, GL_UNSIGNED_BYTE, pixels);

	double modelview_matrix1[16], projection_matrix1[16];
	int viewport1[4];
	glGetDoublev(GL_MODELVIEW_MATRIX, modelview_matrix1);
	glGetDoublev(GL_PROJECTION_MATRIX, projection_matrix1);
	glGetIntegerv(GL_VIEWPORT, viewport1);

	for (int i = 0; i < poly->nquads; i++) { //go through all the quads

		Quad *temp_q = poly->qlist[i];

		glBegin(GL_QUADS);

		for (int j = 0; j < 4; j++) {
			Vertex *temp_v = temp_q->verts[j];

			double x = temp_v->x;
			double y = temp_v->y;

			double tx, ty, dummy;

			gluProject((GLdouble)temp_v->x, (GLdouble)temp_v->y, (GLdouble)temp_v->z,
				modelview_matrix1, projection_matrix1, viewport1, &tx, &ty, &dummy);

			tx = tx / win_width;
			ty = ty / win_height;

			icVector2 dp = icVector2(temp_v->vx, temp_v->vy);
			normalize(dp);

			double dx = dp.x;
			double dy = dp.y;

			double r = dx * dx + dy * dy;
			if (r > dmax*dmax) {
				r = sqrt(r);
				dx *= dmax / r;
				dy *= dmax / r;
			}

			float px = tx + dx;
			float py = ty + dy;

			glTexCoord2f(px, py);
			glVertex3d(temp_v->x, temp_v->y, temp_v->z);
		}
		glEnd();
	}

	iframe = iframe + 1;

	glEnable(GL_BLEND);

	/*blend the drawing with another noise image*/
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();


	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();

	glTranslatef(-1.0, -1.0, 0.0);
	glScalef(2.0, 2.0, 1.0);

	glCallList(iframe % Npat + 1);

	glBegin(GL_QUAD_STRIP);

	glTexCoord2f(0.0, 0.0);  glVertex2f(0.0, 0.0);
	glTexCoord2f(0.0, tmax); glVertex2f(0.0, 1.0);
	glTexCoord2f(tmax, 0.0);  glVertex2f(1.0, 0.0);
	glTexCoord2f(tmax, tmax); glVertex2f(1.0, 1.0);
	glEnd();
	glDisable(GL_BLEND);

	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();

	glMatrixMode(GL_PROJECTION);
	glPopMatrix();

	glReadPixels(0, 0, win_width, win_height, GL_RGB, GL_UNSIGNED_BYTE, pixels);


	/*draw the model with using pixels, note the tx and ty do not take the vector on points*/
	glClearColor(1.0, 1.0, 1.0, 1.0);  // background for rendering color coding and lighting
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, win_width, win_height, 0, GL_RGB, GL_UNSIGNED_BYTE, pixels);
	for (int i = 0; i < poly->nquads; i++) { //go through all the quads
		Quad *temp_q = poly->qlist[i];
		glBegin(GL_QUADS);
		for (int j = 0; j < 4; j++) {
			Vertex *temp_v = temp_q->verts[j];
			double x = temp_v->x;
			double y = temp_v->y;
			double tx, ty, dummy;
			gluProject((GLdouble)temp_v->x, (GLdouble)temp_v->y, (GLdouble)temp_v->z,
				modelview_matrix1, projection_matrix1, viewport1, &tx, &ty, &dummy);
			tx = tx / win_width;
			ty = ty / win_height;
			glTexCoord2f(tx, ty);
			glVertex3d(temp_v->x, temp_v->y, temp_v->z);
		}
		glEnd();
	}

	glDisable(GL_TEXTURE_2D);
	glShadeModel(GL_SMOOTH);
	glDisable(GL_BLEND);
}

/******************************************************************************
Callback function for scene display
******************************************************************************/

void display(void)
{
	glClearColor(1.0, 1.0, 1.0, 1.0);  // background for rendering color coding and lighting

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	set_view(GL_RENDER);
	CHECK_GL_ERROR();

	set_scene(GL_RENDER, poly);
	CHECK_GL_ERROR();

	/*display the mesh*/
	display_polyhedron(poly);
	CHECK_GL_ERROR();

	/*display selected elements*/
	display_selected_vertex(poly);
	CHECK_GL_ERROR();

	display_selected_quad(poly);
	CHECK_GL_ERROR();

	glFlush();
	glutSwapBuffers();
	glFinish();

	CHECK_GL_ERROR();
}


/******************************************************************************
Process a keyboard action.  In particular, exit the program when an
"escape" is pressed in the window.
******************************************************************************/

/*global variable to save polylines*/
PolyLine pentagon;

void keyboard(unsigned char key, int x, int y) {
	int i;

	/* set escape key to exit */
	switch (key) {
	case 27:
		poly->finalize();  // finalize_everything
		exit(0);
		break;

	case '1':
		display_mode = 1;
		glutPostRedisplay();
		break;

	case '2':
		display_mode = 2;
		glutPostRedisplay();
		break;

	case '3':
	{
		display_mode = 3;

		double L = (poly->radius * 2) / 30;
		for (int i = 0; i < poly->nquads; i++) {
			Quad* temp_q = poly->qlist[i];
			for (int j = 0; j < 4; j++) {

				Vertex* temp_v = temp_q->verts[j];

				temp_v->R = int(temp_v->x / L) % 2 == 0 ? 1 : 0;
				temp_v->G = int(temp_v->y / L) % 2 == 0 ? 1 : 0;
				temp_v->B = 0.0;
			}
		}
		glutPostRedisplay();
	}
	break;
	case '4':
		display_mode = 4;
		{
			//examples for dot drawing and polyline drawing

			//create a polylines of a pentagon
			//clear current polylines
			pentagon.clear();
			//there are five vertices of a pentagon
			//the angle of each edge is 2pi/5.0
			double da = 2.0*PI / 5.0;
			for (int i = 0; i < 5; i++) {
				double angle = i * da;
				double cx = cos(angle);
				double cy = sin(angle);

				double n_angle = (i + 1) % 5 * da;
				double nx = cos(n_angle);
				double ny = sin(n_angle);

				LineSegment line(cx, cy, 0, nx, ny, 0);
				pentagon.push_back(line);
			}

		}
		glutPostRedisplay();
		break;

	case '5':
		display_mode = 5;
		//show the IBFV of the field
		break;


	case 't':
	{
		//test case used for putting color in each vertex
		display_mode = 3;
		double min = poly->vlist[0]->scalar;
		double max = poly->vlist[0]->scalar;
		for (int i = 0; i < poly->nverts; i++) {
			Vertex* temp_v = poly->vlist[i];
			double scalar = temp_v->scalar;
			if (scalar < min) {
				min = scalar;
			}
			if (scalar > max) {
				max = scalar;
			}

		}
		//double rgb1[] = { 180, 20, 50 };
		//double rgb2[] = { 50,20,180 };
		double rgb1[] = { 180, 20, 50 };
		double rgb2[] = { 50,150,30 };

		std::vector<int> c1 = { 255,0,255 };
		std::vector<int> c2 = { 153,255,51 };
		for (int i = 0; i < poly->nverts; i++) {
			Vertex* temp_v = poly->vlist[i];
			double scalar =(temp_v->scalar - min) / (max - min);

			double* rgb = interpolateColor(rgb1, rgb2, scalar);

			temp_v->R = rgb[0] / 255;
			temp_v->G = rgb[1] / 255;
			temp_v->B = rgb[2] / 255;

			printf("r: %f\tg: %f\tb: %f\n", rgb[0], rgb[1], rgb[2]);
			
			delete[] rgb;

			//temp_v->R = ((255 / 255 * (scalar - min) / (max - min)) + (153 / 255 * (max - scalar) / (max - min)));
			//temp_v->G = ((0 * (scalar - min) / (max - min)) + (255 / 255 * (max - scalar) / (max - min)));
			//temp_v->B = ((255 / 255 * (scalar - min) / (max - min)) + (51 / 255 * (max - scalar) / (max - min)));
			
		}

		glutPostRedisplay();
		break;

	}

	// test case for displaying color map on mesh with control point stretching the colored section
	case 'y':
	{
		//test case used for putting color in each vertex
		display_mode = 3;
		double min = poly->vlist[0]->scalar;
		double max = poly->vlist[0]->scalar;
		for (int i = 0; i < poly->nverts; i++) {
			Vertex* temp_v = poly->vlist[i];
			double scalar = temp_v->scalar;
			if (scalar < min) {
				min = scalar;
			}
			if (scalar > max) {
				max = scalar;
			}

		}
		//double rgb1[] = { 180, 20, 50 };
		//double rgb2[] = { 50,20,180 };
		double rgb1[] = { 180, 20, 50 };
		double rgb3[] = { 130, 60, 45 }; // control point
		double rgb2[] = { 50,150,30 };

		double control_point_loc = 0.5; // range in [0,1] to put the control point. Scalar values at this point will be the same color as control point
										// values between 0 and this value will be interpolated between rgb1 and rbg3. values above interpolate between rgb3 and rgb2

		std::vector<int> c1 = { 255,0,255 };
		std::vector<int> c2 = { 153,255,51 };
		for (int i = 0; i < poly->nverts; i++) {
			Vertex* temp_v = poly->vlist[i];
			double scalar = (temp_v->scalar - min) / (max - min);
			double* rgb;
			if (scalar <= control_point_loc) {
				scalar = (scalar) / (control_point_loc);
				rgb = interpolateColor(rgb1, rgb3, scalar);
			}
			else {
				scalar = (scalar - control_point_loc) / (1 - control_point_loc);
				rgb = interpolateColor(rgb3, rgb2, scalar);
			}


			temp_v->R = rgb[0] / 255;
			temp_v->G = rgb[1] / 255;
			temp_v->B = rgb[2] / 255;

			printf("r: %f\tg: %f\tb: %f\n", rgb[0], rgb[1], rgb[2]);

			delete[] rgb;

			//temp_v->R = ((255 / 255 * (scalar - min) / (max - min)) + (153 / 255 * (max - scalar) / (max - min)));
			//temp_v->G = ((0 * (scalar - min) / (max - min)) + (255 / 255 * (max - scalar) / (max - min)));
			//temp_v->B = ((255 / 255 * (scalar - min) / (max - min)) + (51 / 255 * (max - scalar) / (max - min)));

		}

		glutPostRedisplay();
		break;

	}

	case 'u':
	{
		get_contours(0);
		display_mode = 6;
		break;
	}
	
	//press 't' to color the space and this case to test
	case'p':
	{
		//case to print the color of the selected point
		showPickedPoint = !showPickedPoint;
		if (showPickedPoint) {
			//get the dimension
			for (int i = 0; i < poly->nverts; i++) {
				if (i == 0) {
					//create minx, minym maxx, maxy variables in polyhedron to save the dimension
					poly->minx = poly->vlist[i]->x;
					poly->maxx = poly->vlist[i]->x;
					poly->miny = poly->vlist[i]->y;
					poly->maxy = poly->vlist[i]->y;
				}
				else {
					if (poly->vlist[i]->x < poly->minx)
						poly->minx = poly->vlist[i]->x;
					if (poly->vlist[i]->x > poly->maxx)
						poly->maxx = poly->vlist[i]->x;
					if (poly->vlist[i]->y < poly->miny)
						poly->miny = poly->vlist[i]->y;
					if (poly->vlist[i]->y > poly->maxy)
						poly->maxy = poly->vlist[i]->y;
				}
			}
			printf("The x coordinate of mesh ranges [%.3f,%.3f]\n", poly->minx, poly->maxx);
			printf("The y coordinate of mesh ranges [%.3f,%.3f]\n", poly->miny, poly->maxy);

			//type the picked point
			float input_x, input_y;
			bool valid_input_x = false;
			bool valid_input_y = false;

			cout << "please input the location" <<endl;
			while (!valid_input_x) {
				cout << "please input the x:" << endl;
				cin >> input_x;
				if (input_x <= poly->maxx && input_x >= poly->minx) {
					valid_input_x = true;
				}
			}
			while (!valid_input_y) {
				cout << "please input the y:" <<endl;
				cin >> input_y;
				if (input_y <= poly->maxy && input_y >= poly->miny) {
					valid_input_y = true;
				}
			}
			pickedPoint.set(input_x, input_y, 0);
			//print_rgb(get_vertex(input_x, input_y));
			
			//test getting the LAB values
			//double* lab = xyz2lab(get_vertex(input_x, input_y));
			//printf("printing LAB values for vertex\n");
			//printf("L: %f, A: %f, B: %f\n", lab[0], lab[1], lab[2]);

			double* xyz_arr1 = rgb2xyz(180, 20, 180);
			printf("hello x: %f, y: %f, z: %f\n", xyz_arr1[0], xyz_arr1[1], xyz_arr1[2]);

			double* lab_arr1 = xyz2lab(xyz_arr1[0], xyz_arr1[1], xyz_arr1[2]);
			printf("hello l: %f, a: %f, b: %f\n", lab_arr1[0], lab_arr1[1], lab_arr1[2]);

			double* xyz_arr2 = lab2xyz(lab_arr1[0], lab_arr1[1], lab_arr1[2]);
			printf("hello x: %f, y: %f, z: %f\n", xyz_arr2[0], xyz_arr2[1], xyz_arr2[2]);

			double* rgb_arr = xyz2rgb(xyz_arr2[0], xyz_arr2[1], xyz_arr2[2]);
			printf("hello r: %f, g: %f, b: %f\n", rgb_arr[0], rgb_arr[1], rgb_arr[2]);
			delete[] xyz_arr1;
			delete[] lab_arr1;
			delete[] xyz_arr2;
			delete[] rgb_arr;


			double* msh_array = rgb2msh(get_vertex(input_x, input_y));
			Vertex* v = get_vertex(input_x, input_y);
			printf("actual r: %f, g: %f, b: %f\n", v->R, v->G, v->B);
			printf("actual M: %f, S: %f, H: %f\n", msh_array[0], msh_array[1], msh_array[2]);
			double* rgb_array = msh2rgb(msh_array[0], msh_array[1], msh_array[2]);
			printf("converted back to r: %f, g: %f, b: %f\n", rgb_array[0], rgb_array[1], rgb_array[2]);
			delete[] msh_array;


			double* linearRgb = linear_rgb(180, 20, 180);
			printf("calc linear rgb r: %f, g: %f, b: %f\n", linearRgb[0], linearRgb[1], linearRgb[2]);
			double* sRgb = sRGB(linearRgb[0], linearRgb[1], linearRgb[2]);
			printf("new sRgb back to r: %f, g: %f, b: %f\n", sRgb[0], sRgb[1], sRgb[2]);

			/*double* xyz_arr = rgb2xyz(12.743, 0.69954, 12.743);
			double x = xyz_arr[0];
			double y = xyz_arr[1];
			double z = xyz_arr[2];
			printf("x: %f, y: %f, z: %f\n", xyz_arr[0], xyz_arr[1], xyz_arr[2]);
			double* rgb_arr = xyz2rgb(x, y, z);
			printf("r: %f, g: %f, b: %f\n", rgb_arr[0], rgb_arr[1], rgb_arr[2]);*/

			//printf("M: %f, S: %f, H: %f\n", *msh_array[0],*msh_array[1], *msh_array[2]);


		}
	
		break;
	}

	case 'r':
		mat_ident(rotmat);
		translation[0] = 0;
		translation[1] = 0;
		zoom = 1.0;
		glutPostRedisplay();
		break;
	}
}



/******************************************************************************
Diaplay the polygon with visualization results
******************************************************************************/


void display_polyhedron(Polyhedron* poly)
{
	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1., 1.);

	glEnable(GL_DEPTH_TEST);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glShadeModel(GL_SMOOTH);
	CHECK_GL_ERROR();

	switch (display_mode) {
	case 1:
	{
		glEnable(GL_LIGHTING);
		glEnable(GL_LIGHT0);
		glEnable(GL_LIGHT1);

		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		GLfloat mat_diffuse[4] = { 1.0, 1.0, 0.0, 0.0 };
		GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
		glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
		glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
		glMaterialf(GL_FRONT, GL_SHININESS, 50.0);

		for (int i = 0; i < poly->nquads; i++) {
			Quad* temp_q = poly->qlist[i];
			glBegin(GL_POLYGON);
			for (int j = 0; j < 4; j++) {
				Vertex* temp_v = temp_q->verts[j];
				glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
		}

		CHECK_GL_ERROR();
	}
	break;

	case 2:
	{
		glDisable(GL_LIGHTING);
		glEnable(GL_LINE_SMOOTH);
		glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glLineWidth(1.0);
		for (int i = 0; i < poly->nquads; i++) {
			Quad* temp_q = poly->qlist[i];

			glBegin(GL_POLYGON);
			for (int j = 0; j < 4; j++) {
				Vertex* temp_v = temp_q->verts[j];
				glNormal3d(temp_q->normal.entry[0], temp_q->normal.entry[1], temp_q->normal.entry[2]);
				glColor3f(0.0, 0.0, 0.0);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();

		}

		glDisable(GL_BLEND);
	}
	break;

	case 3:
		glDisable(GL_LIGHTING);
		for (int i = 0; i < poly->nquads; i++) {
			Quad* temp_q = poly->qlist[i];
			glBegin(GL_POLYGON);
			for (int j = 0; j < 4; j++) {
				Vertex* temp_v = temp_q->verts[j];
				glColor3f(temp_v->R, temp_v->G, temp_v->B);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
		}
		break;

	case 4:
	{
		//draw a dot at position (0.2, 0.3, 0.4) 
		//with radius 0.1 in color blue(0.0, 0.0, 1.0)
		drawDot(0.2, 0.3, 0.4, 0.1, 0.0, 0.0, 1.0);

		//draw a dot at position of vlist[110]
		//with radius 0.2 in color magenta (1.0, 0.0, 1.0)
		Vertex *v = poly->vlist[110];
		drawDot(v->x, v->y, v->z, 0.2, 1.0, 0.0, 1.0);

		//draw line segment start at vlist[110] and end at (vlist[135]->x, vlist[135]->y, 4)
		//with color (0.02, 0.1, 0.02) and width 1
		LineSegment line(poly->vlist[110]->x, poly->vlist[110]->y, poly->vlist[110]->z,
			poly->vlist[135]->x, poly->vlist[135]->y, 4);
		drawLineSegment(line, 1.0, 0.0, 1.0, 0.0);

		//draw a polyline of pentagon with color orange(1.0, 0.5, 0.0) and width 2
		drawPolyline(pentagon, 2.0, 1.0, 0.5, 0.0);

		//display the mesh with color cyan (0.0, 1.0, 1.0)
		glDisable(GL_LIGHTING);
		for (int i = 0; i < poly->nquads; i++) {
			Quad* temp_q = poly->qlist[i];
			glBegin(GL_POLYGON);
			for (int j = 0; j < 4; j++) {
				Vertex* temp_v = temp_q->verts[j];
				glColor3f(0.0, 1.0, 1.0);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
		}
	}
	break;

	case 6:
	{
		double scalar_values[num_intervals] = {};
		for (int i = 0; i < num_intervals; i++) {
			scalar_values[i] = contours[i].second;
		}

		double max = scalar_values[num_intervals - 1];
		double min = scalar_values[0];

		double rgb1[] = { 180, 20, 50 };
		double rgb2[] = { 50,150,30 };

		for (int i = 0; i < num_intervals; i++) {
			double scalar = (scalar_values[i] - min) / (max - min);

			double* rgb = interpolateColor(rgb1, rgb2, scalar);

			drawPolyline((contours[i].first), 2.0, rgb[0] / 255, rgb[1] / 255, rgb[2] / 255);
			delete[] rgb;
		}

		//display the mesh with color cyan (0.0, 1.0, 1.0)
		glDisable(GL_LIGHTING);
		for (int i = 0; i < poly->nquads; i++) {
			Quad* temp_q = poly->qlist[i];
			glBegin(GL_POLYGON);
			for (int j = 0; j < 4; j++) {
				Vertex* temp_v = temp_q->verts[j];
				glColor3f(0.7, 0.7, 0.7);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
		}
		break;
	}


	case 5:
		displayIBFV();
		break;
	}
}
