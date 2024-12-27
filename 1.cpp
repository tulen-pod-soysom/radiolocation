#include <cstdlib>
#include <iostream>
#include <numeric>
#include <ostream>
#include <random>
#include <fstream>

#include "hooke_jeeves.hpp"
#include "NewMHJ.hpp"


constexpr double speed_of_light = 2.99792458e8;

template <typename T, typename Vec>
auto doppler_shift(T frequency, Vec receiver_vel, Vec emitter_vel,
                   Vec receiver_pos, Vec emitter_pos) {
  Vec dpos = receiver_pos;

  for (auto i = 0; i < dpos.size(); ++i) {
    dpos[i] -= emitter_pos[i];
  }
  auto norm = (double)sqrt(
      std::inner_product(dpos.begin(), dpos.end(), dpos.begin(), 0.0));
  for (auto i = 0; i < dpos.size(); ++i) {
    dpos[i] /= norm;
  }

  double v_rec = std::inner_product(receiver_vel.begin(), receiver_vel.end(),
                                    dpos.begin(), 0.0);
  double v_emt = std::inner_product(emitter_vel.begin(), emitter_vel.end(),
                                    dpos.begin(), 0.0);

  return frequency * (1 + v_rec / speed_of_light) /
         (1 - v_emt / speed_of_light);
}

template <typename T> auto sqr(T a) { return a * a; }

auto functional(std::vector<double> pos, std::vector<double> &rec1_pos,
                std::vector<double> &rec2_pos, std::vector<double> &rec3_pos,
                std::vector<double> &rec4_pos, std::vector<double> &rec1_vel,
                std::vector<double> &rec2_vel, std::vector<double> &rec3_vel,
                std::vector<double> &rec4_vel, std::vector<double> dw,
                double w1) {
  auto &x = pos[0];
  auto &y = pos[1];
  auto &z = pos[2];

  auto &x1 = rec1_pos[0];
  auto &y1 = rec1_pos[1];
  auto &z1 = rec1_pos[2];

  auto &x2 = rec2_pos[0];
  auto &y2 = rec2_pos[1];
  auto &z2 = rec2_pos[2];

  auto &x3 = rec3_pos[0];
  auto &y3 = rec3_pos[1];
  auto &z3 = rec3_pos[2];

  auto &x4 = rec4_pos[0];
  auto &y4 = rec4_pos[1];
  auto &z4 = rec4_pos[2];

  auto &vx1 = rec1_vel[0];
  auto &vy1 = rec1_vel[1];
  auto &vz1 = rec1_vel[2];

  auto &vx2 = rec2_vel[0];
  auto &vy2 = rec2_vel[1];
  auto &vz2 = rec2_vel[2];

  auto &vx3 = rec3_vel[0];
  auto &vy3 = rec3_vel[1];
  auto &vz3 = rec3_vel[2];

  auto &vx4 = rec4_vel[0];
  auto &vy4 = rec4_vel[1];
  auto &vz4 = rec4_vel[2];

  auto dw12 = dw[0];
  auto dw13 = dw[1];
  auto dw14 = dw[2];

  auto r1_n = sqrt(sqr(x - x1) + sqr(y - y1) + sqr(z - z1));
  auto r2_n = sqrt(sqr(x - x2) + sqr(y - y2) + sqr(z - z2));
  auto r3_n = sqrt(sqr(x - x3) + sqr(y - y3) + sqr(z - z3));
  auto r4_n = sqrt(sqr(x - x4) + sqr(y - y4) + sqr(z - z4));

  auto v_r1 = (vx1 * (x - x1) + vy1 * (y - y1) + vz1 * (z - z1)) / r1_n;
  auto v_r2 = (vx2 * (x - x2) + vy2 * (y - y2) + vz2 * (z - z2)) / r2_n;
  auto v_r3 = (vx3 * (x - x3) + vy3 * (y - y3) + vz3 * (z - z3)) / r3_n;
  auto v_r4 = (vx4 * (x - x4) + vy4 * (y - y4) + vz4 * (z - z4)) / r4_n;

  auto f1 = (v_r1 - v_r2) / (speed_of_light + v_r1) - dw12 / w1;
  auto f2 = (v_r1 - v_r3) / (speed_of_light + v_r1) - dw13 / w1;
  auto f3 = (v_r1 - v_r4) / (speed_of_light + v_r1) - dw14 / w1;

  return sqr(f1) + sqr(f2) + sqr(f3);
}

struct mat_dot {
  // boost::numeric::ublas::vector<double> pos, vel;
  std::vector<double> pos, vel;

  mat_dot(double x = 0, double y = 0, double z = 0, double vx = 0,
          double vy = 0, double vz = 0)
      : pos(3), vel(3) {
    pos.at(0) = x;
    pos.at(1) = y;
    pos.at(2) = z;
    vel.at(0) = vx;
    vel.at(1) = vy;
    vel.at(2) = vz;
  };
};

std::random_device rd;

auto generate_receiver(double R, double v) {
  double angle =
      std::uniform_real_distribution<double>(0, 0.5 * 3.14159265358979323)(rd);
  double vel_angle =
      std::uniform_real_distribution<double>(0, 2 * 3.14159265358979323)(rd);

  double x = R * cos(angle);
  double y = R * sin(angle);
  double z = 0;
  double vx = v * sqrt(cos(vel_angle) * cos(vel_angle) / (1 + x * x / y / y));
  double vy = v * sqrt(cos(vel_angle) * cos(vel_angle) / (1 + y * y / x / x));
  double vz = v * sin(vel_angle);

  return mat_dot{x, y, z, vx, vy, vz};
}

auto generate_emitter(double R) {
  double angle =
      std::uniform_real_distribution<double>(0, 0.5 * 3.14159265358979323)(rd);
  return mat_dot{R * cos(angle), R * sin(angle), 0};
}

auto print(std::ostream &str, mat_dot &md) {

  str << "pos:" << md.pos.at(0) << '\t' << md.pos.at(1) << '\t' << md.pos.at(2)
      << std::endl;
  str << "vel:" << md.vel.at(0) << '\t' << md.vel.at(1) << '\t' << md.vel.at(2)
      << std::endl;
}

template <typename T>
auto print(std::ostream &str, std::vector<T> &v) {
  for (auto i = 0; i < v.size(); ++i)
  {
    str << v[i] << ' ';
  } str << std::endl;
}

int main() {
  mat_dot receiver[4];
  receiver[0] = generate_receiver(500, 1);
  receiver[1] = generate_receiver(550, 1);
  receiver[2] = generate_receiver(600, 1);
  receiver[3] = generate_receiver(650, 1);

  mat_dot emitter = generate_emitter(100);

  print(std::cout, receiver[0]);
  print(std::cout, receiver[1]);
  print(std::cout, receiver[2]);
  print(std::cout, receiver[3]);

  std::cout << std::endl << "Emitter: " << std::endl;
  print(std::cout, emitter);
  std::cout << std::endl;

  double f = 1e9;

  double f1 = doppler_shift(f, receiver[0].vel, emitter.vel, receiver[0].pos,
                            emitter.pos);
  double f2 = doppler_shift(f, receiver[1].vel, emitter.vel, receiver[1].pos,
                            emitter.pos);
  double f3 = doppler_shift(f, receiver[2].vel, emitter.vel, receiver[2].pos,
                            emitter.pos);
  double f4 = doppler_shift(f, receiver[3].vel, emitter.vel, receiver[3].pos,
                            emitter.pos);

  auto &r1 = receiver[0];
  auto &r2 = receiver[1];
  auto &r3 = receiver[2];
  auto &r4 = receiver[3];

  // auto p = hooke_jeeves_optimizer_simple(
  //     [&](std::vector<double> pos) {
  //       return functional(pos,
  //        r1.pos, r2.pos, r3.pos, r4.pos,
  //        r1.vel, r2.vel,r3.vel, r4.vel,
  //                         std::vector{f2 - f1, f3 - f1, f4 - f1}, f1);
  //     },
  //     /*r1.pos*/std::vector{0.0, 0.0, 0.0}, 3,1e-9);

  std::vector<float> p(3);

  NewMHJ(3,p.data(),[&](float pos[3]) {
               return functional(std::vector<double>(pos,pos+3),
                r1.pos, r2.pos, r3.pos, r4.pos,
                r1.vel, r2.vel,r3.vel, r4.vel,
                                 std::vector{f2 - f1, f3 - f1, f4 - f1}, f1);});


  std::cout << "result:" << p[0] << ' ' << p[1] << ' ' << p[2] << std::endl;


  std::ofstream str("out.txt");
  std::ofstream str2("out_2.txt");

  print(str2,emitter.pos);
  print(str2,p);
  print(str2,r1.pos);
  print(str2,r2.pos);
  print(str2,r3.pos);
  print(str2,r4.pos);

  for (auto j = -1000.; j < 1000.; j += 1)
  {
  for (auto i = -1000.; i < 1000.; i += 1)
  {
    std::vector pos{i,j,0.0};
    str << functional(pos,
                             r1.pos, r2.pos, r3.pos, r4.pos,
                             r1.vel, r2.vel,r3.vel, r4.vel,
                             std::vector{f2 - f1, f3 - f1, f4 - f1}, f1) << ' '/* << i << ' ' << j*/;
  } str << std::endl;
  }
}
