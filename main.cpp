#define _USE_MATH_DEFINES
#include "headers.h"

struct point { double x, y, z; };
struct triangle { vector<point> points; int color; };

inline bool operator==(const point& a, const point& b) {
  return a.x == b.x && a.y == b.y;
}

point get_shape_centroid(const vector<triangle> &shape) {
	double x = 0, y = 0, z = 0;
  int count_side = shape.size();
	for (size_t i = 0; i < count_side; i++) {
		for (size_t j = 0; j < 3; j++) {
			x += shape[i].points[j].x; y += shape[i].points[j].y; z += shape[i].points[j].z;
		}
	}
	return {x /= count_side*3, y /= count_side*3, z /= count_side*3};
}

point get_poly_centroid(const vector<point> &poly) {
	double x = 0, y = 0, z = 0;
  int count = poly.size();
	for (size_t i = 0; i < count; i++) x += poly[i].x, y += poly[i].y, z += poly[i].z;
	return {x /= count, y /= count, z /= count};
}

void rotate(vector<triangle> &shape, double ang[3]) {
	point center = get_shape_centroid(shape);
	Matrix<double> rotate_x(vector<vector<double>>({{1, 0, 0}, {0, cos(ang[0]), sin(ang[0])}, {0, -sin(ang[0]), cos(ang[0])}}));
	Matrix<double> rotate_y(vector<vector<double>>({{cos(ang[1]), 0, -sin(ang[1])}, {0, 1, 0}, {sin(ang[1]), 0, cos(ang[1])}}));
	Matrix<double> rotate_z(vector<vector<double>>({{cos(ang[2]), sin(ang[2]), 0}, {-sin(ang[2]), cos(ang[2]), 0}, {0, 0, 1}}));
	Matrix<double> res = rotate_x*rotate_y*rotate_z;

	double x, y, z;
	for (size_t i = 0; i < shape.size(); i++) {
		for (size_t j = 0; j < 3; j++) {
			x = shape[i].points[j].x; y = shape[i].points[j].y; z = shape[i].points[j].z;
			Matrix<double> temp(vector<vector<double>>({{x-center.x, y-center.y, z-center.z}}));
			Matrix<double> new_coords = temp*res;
			shape[i].points[j].x = new_coords(0, 0) + center.x;
			shape[i].points[j].y = new_coords(0, 1) + center.y;
			shape[i].points[j].z = new_coords(0, 2) + center.z;
		}
	}
	delete ang;
}

void scale(vector<triangle> &shape, double coef) {
	point center = get_shape_centroid(shape);
	coef = coef > 0 ? coef : -1 / coef;
	Matrix<double> scale_mat(vector<vector<double>>({{coef, 0, 0}, {0, coef, 0}, {0, 0, coef}}));

	double x, y, z;
	for (size_t i = 0; i < shape.size(); i++) {
		for (size_t j = 0; j < 3; j++) {
			x = shape[i].points[j].x; y = shape[i].points[j].y; z = shape[i].points[j].z;
			Matrix<double> temp(vector<vector<double>>({{x-center.x, y-center.y, z-center.z}}));
			Matrix<double> new_coords = temp*scale_mat;
			shape[i].points[j].x = new_coords(0, 0) + center.x;
			shape[i].points[j].y = new_coords(0, 1) + center.y;
			shape[i].points[j].z = new_coords(0, 2) + center.z;
		}
	}
}

void translate(vector<triangle> &shape, double dx, double dy, double dz) {
	for (size_t i = 0; i < shape.size(); i++) {
		for (size_t j = 0; j < 3; j++) {
			shape[i].points[j].x += dx, shape[i].points[j].y += dy, shape[i].points[j].z += dz;
		}
	}
}

bool is_intersection(const vector<point> line_p, const vector<point> line_w, double *x, double *y) {
  double c1_x, c1_y, c2_x, c2_y, s, t;
	c1_x = line_p[1].x - line_p[0].x;
	c1_y = line_p[1].y - line_p[0].y;
	c2_x = line_w[1].x - line_w[0].x;
	c2_y = line_w[1].y - line_w[0].y;

	s = (-c1_y * (line_p[0].x - line_w[0].x) + c1_x * (line_p[0].y - line_w[0].y)) / (-c2_x * c1_y + c1_x * c2_y);
	t = ( c2_x * (line_p[0].y - line_w[0].y) - c2_y * (line_p[0].x - line_w[0].x)) / (-c2_x * c1_y + c1_x * c2_y);
	
	if (s >= 0 && s <= 1 && t >= 0 && t <= 1) {
		*x = line_p[0].x + (t * c1_x); *y = line_p[0].y + (t * c1_y);
		return true;
	}
	return false;
}

bool point_in_window(point p, vector<point> w) {
  return p.x >= w[0].x && p.x <= w[1].x && p.y <= w[2].y && p.y >= w[0].y;
}

double sign (point p1, point p2, point p3) {
  return (p1.x - p3.x) * (p2.y - p3.y) - (p2.x - p3.x) * (p1.y - p3.y);
}

bool point_in_triangle(point p, vector<point> t) {
  double d1, d2, d3;
  d1 = sign(p, t[0], t[1]);
  d2 = sign(p, t[1], t[2]);
  d3 = sign(p, t[2], t[0]);
  return !(((d1 <= 0) || (d2 <= 0) || (d3 <= 0)) && ((d1 >= 0) || (d2 >= 0) || (d3 >= 0)));
}

bool point_in_poly(point p, vector<point> poly) {
  for(point &_p : poly) if (_p == p) return true;
  return false;
}

vector<point> get_polygon_intersection(const vector<point> &p, const vector<point> &w) {
  vector<point> res;
  int count_vert_p = p.size(), count_vert_w = w.size();
  double x, y;

  for (size_t i = 0; i < count_vert_w; i++) {
    if (point_in_triangle(w[i], p)) res.push_back(w[i]);
  }

  for (size_t i = 0; i < count_vert_p; i++) {
    if (point_in_window(p[i], w)) res.push_back(p[i]);
  }

  for (size_t i = 0; i < count_vert_p; i++) {
    vector<point> line_p = {p[i], p[(i+1) % count_vert_p]};
    for (size_t j = 0; j < count_vert_w; j++) {
      vector<point> line_w = {w[j], w[(j+1) % count_vert_w]};
      
      if (is_intersection(line_p, line_w, &x, &y)) {
        if (!point_in_poly({x, y}, res)) res.push_back({x, y});
      } 
    }
  }

  point center = get_poly_centroid(res);
  sort(res.begin(), res.end(), [center](point a, point b) {
    if (a.x - center.x >= 0 && b.x - center.x < 0) return true;
    if (a.x - center.x < 0 && b.x - center.x >= 0) return false;
    if (a.x - center.x == 0 && b.x - center.x == 0) {
      if (a.y - center.y >= 0 || b.y - center.y >= 0) return a.y > b.y;
      return b.y > a.y;
    }

    int det = (a.x - center.x) * (b.y - center.y) - (b.x - center.x) * (a.y - center.y);
    if (det < 0) return true;
    if (det > 0) return false;

    int d1 = (a.x - center.x) * (a.x - center.x) + (a.y - center.y) * (a.y - center.y);
    int d2 = (b.x - center.x) * (b.x - center.x) + (b.y - center.y) * (b.y - center.y);
    return d1 > d2;
  });
  return res;
}

void fill(const vector<point> &poly, int color=WHITE) {
	setcolor(color);
	int count_vertex = poly.size();
	double bounds[2][2] = {{poly[0].x, poly[0].y}, {poly[0].x, poly[0].y}};
	for (int i = 1; i < count_vertex; i++) {
		int x = poly[i].x, y = poly[i].y;
		if (x < bounds[0][0]) bounds[0][0] = x;
		else if (x > bounds[1][0]) bounds[1][0] = x;
		if (y < bounds[0][1]) bounds[0][1] = y;
		else if (y > bounds[1][1]) bounds[1][1] = y;
	}
	for (double i = bounds[0][1]+ 0.01; i <= bounds[1][1]; i++) {
		double pts[20] = {0};
		int count = 0;

		vector<point> pt_vec = {{bounds[0][0] - 1, i}, {bounds[1][0] + 1, i}};
		for (int j = 0; j < count_vertex; j++)
		{
			vector<point> line = {{poly[j].x, poly[j].y},
													{poly[(j + 1) % count_vertex].x, poly[(j + 1) % count_vertex].y}};
			double x, y;
			if (is_intersection(pt_vec, line, &x, &y)) pts[count++] = x;
		}
		if (count) {
			qsort(pts, count, sizeof(double), [](const void* x, const void* y) { return (int)(*(double*)x - *(double*)y); });
			for (int j = 0; j < count; j++) {
				if (pts[j+1]) line(pts[j], i, pts[j + 1], i);
			}
		}
	}
}

void show_poly(const vector<point> &poly, int fill_color=WHITE, int border_color=WHITE) {
  int count_vert = poly.size();
	if (!count_vert) return;
	if (fill_color != -1) fill(poly, fill_color);
	if (border_color != -1) {
		setcolor(border_color);
		for (size_t i = 0; i < count_vert; i++) {
			line(poly[i].x, poly[i].y, poly[(i+1) % count_vert].x, poly[(i+1) % count_vert].y);
		}
	}
}

int get_count_poly_in_window(const vector<triangle> &shape, const vector<point> &win, vector<triangle> &ps) {
	int count = 0;
	for (size_t i = 0; i < shape.size(); i++) {
		vector<point> temp = get_polygon_intersection(shape[i].points, win);
    if (temp.size()) {
			count++;
			ps.push_back({temp, shape[i].color});
		}
	}
	return count;
}

bool is_poly_in_window(vector<point> &poly, vector<point> &win) {
	for (size_t i = 0; i < poly.size(); i++) {
		if (!point_in_window(poly[i], win)) return false;
	}
	return true;
}

void divide(const vector<triangle> &shape, double x1=0, double y1=0, double x2=width, double y2=height) {
	if (x2 - x1 < eps || y2 - y1 < eps) return;
	
	vector<point> window = {{x1, y1}, {x2, y1}, {x2, y2}, {x1, y2}};
  // rectangle(x1, y1, x2, y2);
	vector<triangle> polygons;
	int count = get_count_poly_in_window(shape, window, polygons);

	if (!count) return;
	else if (count == 1) {
		show_poly(polygons[0].points, polygons[0].color, -1);
		return;
	}
	else {
		bool is_window_in_poly;
		vector<triangle> out_tri;

		for (size_t i = 0; i < shape.size(); i++) {
			is_window_in_poly = true;
			for (size_t j = 0; j < 4; j++) {
				if (!point_in_triangle(window[j], shape[i].points)) is_window_in_poly = false;
			}
			if (is_window_in_poly) {
				out_tri.push_back(shape[i]);
			}
		}

		if(out_tri.size() == count) {
			triangle max_z_poly = out_tri[0];
			double avg_max_z = -9999.0, sum_max = 0;

			for (size_t i = 0; i < out_tri.size(); i++) {
				for (size_t j = 0; j < 3; j++) {
					sum_max += out_tri[i].points[j].z;
				}
				if (avg_max_z < sum_max/3.0) {
					avg_max_z = sum_max/3.0;
					max_z_poly = out_tri[i];
				}
				sum_max = 0;
			}
			
			window = {{x1-1, y1-1}, {x2+1, y1-1}, {x2+1, y2+1}, {x1-1, y2+1}};
			show_poly(window, max_z_poly.color, max_z_poly.color);
			return;
		}
	}
  double xc = (x2 - x1) / 2 + x1;
  double yc = (y2 - y1) / 2 + y1;

  divide(shape, x1, y1, xc, yc);
  divide(shape, x1, yc, xc, y2);
  divide(shape, xc, y1, x2, yc);
  divide(shape, xc, yc, x2, y2);
}

void show(const vector<vector<triangle>> &shapes) {
	vector<triangle> res;
	for (auto &shape : shapes) res.insert(res.end(), shape.begin(), shape.end());

	int count_side = res.size();
  triangle poly;
	poly.points.resize(3);
	vector<triangle> vec_tri;

  for (size_t i = 0; i < count_side; i++) {
    for (size_t j = 0; j < 3; j++) {
      poly.points[j].x = (res[i].points[j].x * 0.5 + 0.5) * width;
      poly.points[j].y = (res[i].points[j].y * 0.5 + 0.5) * height;
			poly.points[j].z = res[i].points[j].z;
    }
		poly.color = res[i].color;
		vec_tri.push_back(poly);
  }
	
	double x1 = 0, y1 = 0, x2 = width, y2 = height;
  vector<point> window = {{x1, y1}, {x2, y1}, {x2, y2}, {x1, y2}};
	divide(vec_tri, x1,y1,x2,y2);
}

void click_handler(int x, int y) { cout << x << " " << y << "\n"; }

int main() {
	int win = initwindow(width, height, "cg"), key;

	vector<triangle> prism1 = {
		{ {{-0.5, -0.5, -2.0}, {0.0, -0.5, -2.8}, {0.5, -0.5, -2.0}}, 1}, // abc
		{ {{-0.5,  0.5, -2.0}, {0.5,  0.5, -2.0}, {0.0,  0.5, -2.8}}, 2}, // def
		{ {{-0.5, -0.5, -2.0}, {0.5, -0.5, -2.0}, {0.5,  0.5, -2.0}}, 4}, // ace
		{ {{-0.5, -0.5, -2.0}, {-0.5,  0.5, -2.0}, {0.5,  0.5, -2.0}}, 4}, // ade
		{ {{-0.5, -0.5, -2.0}, {0.0, -0.5, -2.8}, {-0.5,  0.5, -2.0}}, 3}, // abd
		{ {{0.0, -0.5, -2.8}, {0.0,  0.5, -2.8}, {-0.5,  0.5, -2.0}}, 3}, // bfd
		{ {{0.0, -0.5, -2.8}, {0.0,  0.5, -2.8}, {0.5,  0.5, -2.0}}, 5}, // bef
		{ {{0.0, -0.5, -2.8}, {0.5, -0.5, -2.0}, {0.5,  0.5, -2.0}}, 5}, // bce
	};

	vector<triangle> prism2 = {
		{ {{-0.5, -0.5, -2.0}, {0.0, -0.5, -2.8}, {0.5, -0.5, -2.0}}, 1}, // abc
		{ {{-0.5,  0.5, -2.0}, {0.5,  0.5, -2.0}, {0.0,  0.5, -2.8}}, 2}, // def
		{ {{-0.5, -0.5, -2.0}, {0.5, -0.5, -2.0}, {0.5,  0.5, -2.0}}, 4}, // ace
		{ {{-0.5, -0.5, -2.0}, {-0.5,  0.5, -2.0}, {0.5,  0.5, -2.0}}, 4}, // ade
		{ {{-0.5, -0.5, -2.0}, {0.0, -0.5, -2.8}, {-0.5,  0.5, -2.0}}, 3}, // abd
		{ {{0.0, -0.5, -2.8}, {0.0,  0.5, -2.8}, {-0.5,  0.5, -2.0}}, 3}, // bfd
		{ {{0.0, -0.5, -2.8}, {0.0,  0.5, -2.8}, {0.5,  0.5, -2.0}}, 5}, // bef
		{ {{0.0, -0.5, -2.8}, {0.5, -0.5, -2.0}, {0.5,  0.5, -2.0}}, 5}, // bce
	};

	scale(prism1, -inc_coef*17);
	scale(prism2, -inc_coef*17);

	translate(prism1, speed*5, speed*5, 0);
	rotate(prism1, new double[3]{deg2rad(40), deg2rad(60), deg2rad(0)});
	rotate(prism2, new double[3]{deg2rad(-40), deg2rad(-60), deg2rad(0)});

	vector<vector<triangle>> shapes = {prism1, prism2};
  setbkcolor(background_color);
	show(shapes);

	bool is_1 = true;
	while (1) {
		key = getch();
		if (key == S_KEY) rotate(is_1 ? shapes[0] : shapes[1], new double[3]{deg2rad(angle), deg2rad(0), deg2rad(0)});
		else if (key == W_KEY) rotate(is_1 ? shapes[0] : shapes[1], new double[3]{deg2rad(-angle), deg2rad(0), deg2rad(0)});
		else if (key == D_KEY) rotate(is_1 ? shapes[0] : shapes[1], new double[3]{deg2rad(0), deg2rad(angle), deg2rad(0)});
		else if (key == A_KEY) rotate(is_1 ? shapes[0] : shapes[1], new double[3]{deg2rad(0), deg2rad(-angle), deg2rad(0)});
		else if (key == X_KEY) rotate(is_1 ? shapes[0] : shapes[1], new double[3]{deg2rad(0), deg2rad(0), deg2rad(angle)});
		else if (key == Z_KEY) rotate(is_1 ? shapes[0] : shapes[1], new double[3]{deg2rad(0), deg2rad(0), deg2rad(-angle)});
		else if (key == PLUS_KEY) scale(is_1 ? shapes[0] : shapes[1], inc_coef);
		else if (key == MINUS_KEY) scale(is_1 ? shapes[0] : shapes[1], -inc_coef);
		else if (key == UP_KEY) translate(is_1 ? shapes[0] : shapes[1], 0, -speed, 0);
		else if (key == DOWN_KEY) translate(is_1 ? shapes[0] : shapes[1], 0, speed, 0);
		else if (key == LEFT_KEY) translate(is_1 ? shapes[0] : shapes[1], -speed, 0, 0);
		else if (key == RIGHT_KEY) translate(is_1 ? shapes[0] : shapes[1], speed, 0, 0);
		else if (key == 91) translate(is_1 ? shapes[0] : shapes[1], 0, 0, -speed);
		else if (key == 93) translate(is_1 ? shapes[0] : shapes[1], 0, 0, speed);
		else if (key == ESC_KEY) is_1 = !is_1;

		show(shapes); swapbuffers(); cleardevice();
	}
	return 0;
}