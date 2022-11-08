#include <iostream>
#include "graphics.h"
#include <vector>
#include <set>
#include <algorithm>
#include <cmath>
#include <array>

using std::vector;
using std::atan2;

#define width 512
#define height 512
#define eps 30

int color = 1;

struct point3d { float x, y, z, angle; };
struct triangle { point3d tri[3]; int color; };
struct shape { vector<triangle> triangles; };

inline bool operator==(const point3d& a, const point3d& b) {
  return a.x == b.x && a.y == b.y;
}

enum location_poly {in, out, intersect, cover};

vector<point3d> triangle_coords = {{100, 100}, {200, 100}, {150, 200}};

bool point_in_window(point3d p, vector<point3d> w) {
  return p.x >= w[0].x && p.x <= w[1].x && p.y <= w[2].y && p.y >= w[0].y;
}

float sign (point3d p1, point3d p2, point3d p3) {
  return (p1.x - p3.x) * (p2.y - p3.y) - (p2.x - p3.x) * (p1.y - p3.y);
}

bool point_in_triangle(point3d p, vector<point3d> t) {
  float d1, d2, d3;
  d1 = sign(p, t[0], t[1]);
  d2 = sign(p, t[1], t[2]);
  d3 = sign(p, t[2], t[0]);
  return !(((d1 < 0) || (d2 < 0) || (d3 < 0)) && ((d1 > 0) || (d2 > 0) || (d3 > 0)));
}

bool point_in_poly(point3d point, vector<point3d> poly) {
  for(point3d &p : poly) if (p == point) return true;
  return false;
}

point3d get_centroid(const vector<point3d> &poly) {
	float x = 0, y = 0, z = 0;
  int count = poly.size();
	for (size_t i = 0; i < count; i++) x += poly[i].x, y += poly[i].y, z += poly[i].z;
	return {x /= count, y /= count, z /= count};
}

vector<point3d> intersection(vector<point3d> &p, const vector<point3d> w) {
  vector<point3d> res;
  int count_vert_p = p.size(), count_vert_w = w.size();
  float c1_x, c1_y, c2_x, c2_y, s, t, x, y;

  for (size_t i = 0; i < count_vert_w; i++) {
    if (point_in_triangle(w[i], p)) res.push_back(w[i]);
  }

  for (size_t i = 0; i < count_vert_p; i++) {
    c1_x = p[(i+1) % count_vert_p].x - p[i].x;
    c1_y = p[(i+1) % count_vert_p].y - p[i].y;
    for (size_t j = 0; j < count_vert_w; j++) {
      c2_x = w[(j+1) % count_vert_w].x - w[j].x;
      c2_y = w[(j+1) % count_vert_w].y - w[j].y;
      s = (-c1_y * (p[i].x - w[j].x) + c1_x * (p[i].y - w[j].y)) / (-c2_x * c1_y + c1_x * c2_y);
	    t = ( c2_x * (p[i].y - w[j].y) - c2_y * (p[i].x - w[j].x)) / (-c2_x * c1_y + c1_x * c2_y);
      if (s >= 0 && s <= 1 && t >= 0 && t <= 1) {
        x = p[i].x + (t * c1_x);
        y = p[i].y + (t * c1_y);
        if (!point_in_poly({x, y}, res)) res.push_back({x, y});
      }
    }
  }

  for (size_t i = 0; i < count_vert_p; i++) {
    if (point_in_window(p[i], w)) res.push_back(p[i]);
  }

  point3d center = get_centroid(res);
  sort(res.begin(), res.end(), [center](point3d a, point3d b) {
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

bool is_in(vector<point3d> p, vector<point3d> w) {
  int count_vert = p.size();
  for (size_t i = 0; i < count_vert; i++) {
    if (p[i].x > w[1].x || p[i].x < w[0].x || p[i].y > w[2].y || p[i].y < w[0].y) return false;
  }
  return true;
}

bool is_out(vector<point3d> p, vector<point3d> w) {
  return false;
}

location_poly get_location(vector<point3d> p, vector<point3d> w) {
  return location_poly::in;
}

void show_poly(const vector<point3d> &poly) {
  int count_vert = poly.size();
  for (size_t i = 0; i < count_vert; i++) {
    line(poly[i].x, poly[i].y, poly[(i+1) % count_vert].x, poly[(i+1) % count_vert].y);
  }
}

void divide(float x1, float y1, float x2, float y2) {
  if (x2 - x1 <= eps || y2 - y1 <= eps) return;
  setcolor(WHITE);
  rectangle(x1, y1, x2, y2);
  vector<point3d> window = {{x1, y1}, {x2, y1}, {x2, y2}, {x1, y2}};
  setcolor(color++%5+1);
  show_poly(intersection(triangle_coords, window));
  Sleep(500);
  float xc = (x2 - x1) / 2 + x1;
  float yc = (y2 - y1) / 2 + y1;

  divide(x1, y1, xc, yc);
  divide(x1, yc, xc, y2);
  divide(xc, y1, x2, yc);
  divide(xc, yc, x2, y2);
}

void click_handler(int x, int y) {
  std::cout << x << " " << y << " <- point\n";
}

int main() {
  int win = initwindow(width, height, "ddd"), key;
  registermousehandler(WM_LBUTTONDOWN, click_handler);
  // divide(0,0,width,height);
  float x1 = 100;
  float y1 = 120;
  float x2 = 150;
  float y2 = 150;
  vector<point3d> window = {{x1, y1}, {x2, y1}, {x2, y2}, {x1, y2}};
  // setcolor(color++%5+1);
  show_poly(window);
  show_poly(triangle_coords);
  setcolor(YELLOW);
  show_poly(intersection(triangle_coords, window));
  while (1) {
		key = getch();
		if (key == 75) x1-=5, x2-=5;
    else if (key == 77) x1+=5, x2+=5;
    if (key == 72) y1-=5, y2-=5;
    else if (key == 80) y1+=5, y2+=5;
    window = {{x1, y1}, {x2, y1}, {x2, y2}, {x1, y2}};
    cleardevice();
    setlinestyle(SOLID_LINE, 1, 1);
    setcolor(WHITE); show_poly(window); show_poly(triangle_coords);
    setlinestyle(SOLID_LINE, 1, 3);
    vector<point3d> intersect_poly = intersection(triangle_coords, window);
    setcolor(YELLOW); show_poly(intersect_poly);
    for (size_t i = 0; i < intersect_poly.size(); i++) {
      std::cout << intersect_poly[i].x << " " << intersect_poly[i].y << "\n";
    }
    std::cout << "\n\n";
  }
  getch();
  return 0;
}