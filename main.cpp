#define _USE_MATH_DEFINES
#include "headers.h"

struct point { float x, y, z; };
struct triangle { vector<point> points; int color; };

inline bool operator==(const point& a, const point& b) {
  return a.x == b.x && a.y == b.y;
}

point get_shape_centroid(const vector<triangle> &shape) {
	float x = 0, y = 0, z = 0;
  int count_side = shape.size();
	for (size_t i = 0; i < count_side; i++) {
		for (size_t j = 0; j < 3; j++) {
			x += shape[i].points[j].x; y += shape[i].points[j].y; z += shape[i].points[j].z;
		}
	}
	return {x /= count_side*3, y /= count_side*3, z /= count_side*3};
}

point get_poly_centroid(const vector<point> &poly) {
	float x = 0, y = 0, z = 0;
  int count = poly.size();
	for (size_t i = 0; i < count; i++) x += poly[i].x, y += poly[i].y, z += poly[i].z;
	return {x /= count, y /= count, z /= count};
}

void rotate(vector<triangle> &shape, float ang[3]) {
	point center = get_shape_centroid(shape);
	Matrix<float> rotate_x(vector<vector<float>>({{1, 0, 0}, {0, cos(ang[0]), sin(ang[0])}, {0, -sin(ang[0]), cos(ang[0])}}));
	Matrix<float> rotate_y(vector<vector<float>>({{cos(ang[1]), 0, -sin(ang[1])}, {0, 1, 0}, {sin(ang[1]), 0, cos(ang[1])}}));
	Matrix<float> rotate_z(vector<vector<float>>({{cos(ang[2]), sin(ang[2]), 0}, {-sin(ang[2]), cos(ang[2]), 0}, {0, 0, 1}}));
	Matrix<float> res = rotate_x*rotate_y*rotate_z;

	float x, y, z;
	for (size_t i = 0; i < shape.size(); i++) {
		for (size_t j = 0; j < 3; j++) {
			x = shape[i].points[j].x; y = shape[i].points[j].y; z = shape[i].points[j].z;
			Matrix<float> temp(vector<vector<float>>({{x-center.x, y-center.y, z-center.z}}));
			Matrix<float> new_coords = temp*res;
			shape[i].points[j].x = new_coords(0, 0) + center.x;
			shape[i].points[j].y = new_coords(0, 1) + center.y;
			shape[i].points[j].z = new_coords(0, 2) + center.z;
		}
	}
	delete ang;
}

void scale(vector<triangle> &shape, float coef) {
	point center = get_shape_centroid(shape);
	coef = coef > 0 ? coef : -1 / coef;
	Matrix<float> scale_mat(vector<vector<float>>({{coef, 0, 0}, {0, coef, 0}, {0, 0, coef}}));

	float x, y, z;
	for (size_t i = 0; i < shape.size(); i++) {
		for (size_t j = 0; j < 3; j++) {
			x = shape[i].points[j].x; y = shape[i].points[j].y; z = shape[i].points[j].z;
			Matrix<float> temp(vector<vector<float>>({{x-center.x, y-center.y, z-center.z}}));
			Matrix<float> new_coords = temp*scale_mat;
			shape[i].points[j].x = new_coords(0, 0) + center.x;
			shape[i].points[j].y = new_coords(0, 1) + center.y;
			shape[i].points[j].z = new_coords(0, 2) + center.z;
		}
	}
}

void translate(vector<triangle> &shape, float dx, float dy, float dz) {
	for (size_t i = 0; i < shape.size(); i++) {
		for (size_t j = 0; j < 3; j++) {
			shape[i].points[j].x += dx, shape[i].points[j].y += dy, shape[i].points[j].z += dz;
		}
	}
}

bool is_intersection(const vector<point> line_p, const vector<point> line_w, float *x, float *y) {
  float c1_x, c1_y, c2_x, c2_y, s, t;
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

float sign (point p1, point p2, point p3) {
  return (p1.x - p3.x) * (p2.y - p3.y) - (p2.x - p3.x) * (p1.y - p3.y);
}

bool point_in_triangle(point p, vector<point> t) {
  float d1, d2, d3;
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
  float x, y;

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
	float bounds[2][2] = {{poly[0].x, poly[0].y}, {poly[0].x, poly[0].y}};
	for (int i = 1; i < count_vertex; i++) {
		int x = poly[i].x, y = poly[i].y;
		if (x < bounds[0][0]) bounds[0][0] = x;
		else if (x > bounds[1][0]) bounds[1][0] = x;
		if (y < bounds[0][1]) bounds[0][1] = y;
		else if (y > bounds[1][1]) bounds[1][1] = y;
	}
	for (float i = bounds[0][1]+ 0.01; i <= bounds[1][1]; i++) {
		float pts[20] = {0};
		int count = 0;

		vector<point> pt_vec = {{bounds[0][0] - 1, i}, {bounds[1][0] + 1, i}};
		for (int j = 0; j < count_vertex; j++)
		{
			vector<point> line = {{poly[j].x, poly[j].y},
													{poly[(j + 1) % count_vertex].x, poly[(j + 1) % count_vertex].y}};
			float x, y;
			if (is_intersection(pt_vec, line, &x, &y)) pts[count++] = x;
		}
		if (count) {
			qsort(pts, count, sizeof(float), [](const void* x, const void* y) { return (int)(*(float*)x - *(float*)y); });
			for (int j = 0; j < count; j++) {
				if (pts[j+1]) line(pts[j], i, pts[j + 1], i);
			}
		}
	}
}

void show_poly(const vector<point> &poly, bool is_filling, int color=WHITE, int border_color=WHITE) {
  int count_vert = poly.size();
	if (!count_vert) return;
	if (is_filling) fill(poly, color);
	if (border_color == -1) return;
	setcolor(border_color);
  for (size_t i = 0; i < count_vert; i++) {
    line(poly[i].x, poly[i].y, poly[(i+1) % count_vert].x, poly[(i+1) % count_vert].y);
  }
	setcolor(WHITE);
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
	// cout << count << "\n";
	return count;
}

bool is_poly_in_window(vector<point> &poly, vector<point> &win) {
	for (size_t i = 0; i < poly.size(); i++) {
		if (!point_in_window(poly[i], win)) return false;
	}
	return true;
}

void divide(const vector<triangle> &shape, float x1=0, float y1=0, float x2=width, float y2=height) {
  if (x2 - x1 <= eps || y2 - y1 <= eps) return;
	
	vector<point> window = {{x1, y1}, {x2, y1}, {x2, y2}, {x1, y2}};
  // rectangle(x1, y1, x2, y2);
	// show_poly(window, false);
	vector<triangle> polygons;
	int count = get_count_poly_in_window(shape, window, polygons); // poly.size()

	if (!count) return;
	else if (count == 1) { // && is_poly_in_window(polygons[0].points, window)
		show_poly(polygons[0].points, true, polygons[0].color);
		return;
	}
	else {
		bool in_window; // охватывает окно
		vector<triangle> out_tri;

		for (size_t i = 0; i < shape.size(); i++) {
			in_window = true;
			for (size_t j = 0; j < 4; j++) {
				if (!point_in_triangle(window[j], shape[i].points)) in_window = false;
				// if (window[j].x == )
			}
			if (in_window) {
				out_tri.push_back(shape[i]);
			}
		}
		// cout << out_tri.size() << " size\n";
		if(out_tri.size()) {
			triangle max_z_poly = out_tri[0];
			float max_z = max_z_poly.points[0].z;

			for (size_t i = 0; i < out_tri.size(); i++) {
				for (size_t j = 0; j < 3; j++) {
					if (out_tri[i].points[j].z > max_z) {
						max_z = out_tri[i].points[j].z;
						max_z_poly = out_tri[i];
					}
				}
			}
			show_poly(window, true, max_z_poly.color, max_z_poly.color);
			// setfillstyle(SOLID_FILL, max_z_poly.color);
			// setcolor(max_z_poly.color);
			// int poly[10] = {x1,y1,x2,y1,x2,y2,x1,y2,x1,y1};
			// drawpoly(5, poly);
			// fillpoly(5, poly);
			return;
		}
	}
	// Sleep(500);
  float xc = (x2 - x1) / 2 + x1;
  float yc = (y2 - y1) / 2 + y1;

  divide(shape, x1, y1, xc, yc);
  divide(shape, x1, yc, xc, y2);
  divide(shape, xc, y1, x2, yc);
  divide(shape, xc, yc, x2, y2);
}

void show(const vector<triangle> &shape) {
	int count_side = shape.size();
  triangle poly; // ???
	poly.points.resize(3);
	vector<triangle> vec_tri;

  for (size_t i = 0; i < count_side; i++) {
    for (size_t j = 0; j < 3; j++) {
      poly.points[j].x = (shape[i].points[j].x * 0.5 + 0.5) * width;
      poly.points[j].y = (shape[i].points[j].y * 0.5 + 0.5) * height;
    }
		poly.color = shape[i].color;
		vec_tri.push_back(poly);
		show_poly(poly.points, true, poly.color);
		// float x1 = 50, y1 = 50, x2 = 390, y2 = 390;
  	// vector<point> window = {{x1, y1}, {x2, y1}, {x2, y2}, {x1, y2}};
		// show_poly(window, false);
		// cout << get_count_poly_in_window(vec_tri, window) << "\n";
		// show_poly(get_polygon_intersection(poly.points, window), true, poly.color);
  }
	float x1 = 250, y1 = 250, x2 = 300, y2 = 300;
  vector<point> window = {{x1, y1}, {x2, y1}, {x2, y2}, {x1, y2}};
	divide(vec_tri, x1,y1,x2,y2);
	show_poly(window, false);
}

int main() {
	int win = initwindow(width, height, "cg"), key;
	// vector<point> window = {{x1, y1}, {x2, y1}, {x2, y2}, {x1, y2}};
	
	vector<triangle> prism1 = {
		{ {{-0.5, -0.5, -2.0}, {0.0, -0.5, -2.8}, {0.5, -0.5, -2.0}}, 1}, // abc
		{ {{-0.5,  0.5, -2.0}, {0.5,  0.5, -2.0}, {0.0,  0.5, -2.8}}, 2}, // def
		{ {{-0.5, -0.5, -2.0}, {0.5, -0.5, -2.0}, {0.5,  0.5, -2.0}}, 3}, // ace
		{ {{-0.5, -0.5, -2.0}, {-0.5,  0.5, -2.0}, {0.5,  0.5, -2.0}}, 4}, // ade
		{ {{-0.5, -0.5, -2.0}, {0.0, -0.5, -2.8}, {-0.5,  0.5, -2.0}}, 5}, // abd
		{ {{0.0, -0.5, -2.8}, {0.0,  0.5, -2.8}, {-0.5,  0.5, -2.0}}, 6}, // bfd
		{ {{0.0, -0.5, -2.8}, {0.5, -0.5, -2.0}, {0.5,  0.5, -2.0}}, 7}, // bce
		{ {{0.0, -0.5, -2.8}, {0.5,  0.5, -2.0}, {0.0,  0.5, -2.8}}, 8}, // bef
  };
	
	// scale(prism1, -inc_coef*5);
	// translate(prism1, speed*5, speed*5, 0);

	show(prism1);
	
  setbkcolor(background_color);

	while (1) {
		key = getch();
		if (key == S_KEY) rotate(prism1, new float[3]{deg2rad(angle), deg2rad(0), deg2rad(0)});
		else if (key == W_KEY) rotate(prism1, new float[3]{deg2rad(-angle), deg2rad(0), deg2rad(0)});
		else if (key == D_KEY) rotate(prism1, new float[3]{deg2rad(0), deg2rad(angle), deg2rad(0)});
		else if (key == A_KEY) rotate(prism1, new float[3]{deg2rad(0), deg2rad(-angle), deg2rad(0)});
		else if (key == X_KEY) rotate(prism1, new float[3]{deg2rad(0), deg2rad(0), deg2rad(angle)});
		else if (key == Z_KEY) rotate(prism1, new float[3]{deg2rad(0), deg2rad(0), deg2rad(-angle)});
		else if (key == PLUS_KEY) scale(prism1, inc_coef);
		else if (key == MINUS_KEY) scale(prism1, -inc_coef);
		else if (key == UP_KEY) translate(prism1, 0, -speed, 0);
		else if (key == DOWN_KEY) translate(prism1, 0, speed, 0);
		else if (key == LEFT_KEY) translate(prism1, -speed, 0, 0);
		else if (key == RIGHT_KEY) translate(prism1, speed, 0, 0);
		show(prism1); swapbuffers(); cleardevice();
		Sleep(0);
	}
	return 0;
}