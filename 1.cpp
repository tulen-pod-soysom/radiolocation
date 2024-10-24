#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_expression.hpp>
#include <cstdlib>
#include <iostream>
#include <ostream>
#include <random>

using namespace boost::numeric::ublas;

constexpr double speed_of_light = 2.99792458e8;

template <typename T, typename Vec>
auto doppler_shift(T frequency, Vec receiver_vel, Vec emitter_vel, Vec receiver_pos, Vec emitter_pos)
{
    Vec dpos = receiver_pos - emitter_pos;
    dpos /= (double)sqrt(inner_prod(dpos, dpos));

    double v_rec = inner_prod(receiver_vel,dpos);
    double v_emt = inner_prod(emitter_vel,dpos);

    return frequency * (1 + v_rec / speed_of_light) / (1 - v_emt / speed_of_light);
}

template <typename T>
auto sqr(T a) {return a*a;}


template <typename Real = double>
auto F1(Real dw12, Real w1,
 Real x, Real y, Real z,
 Real x1,Real y1, Real z1, 
 Real x2,Real y2, Real z2, 
 Real vx1, Real vy1, Real vz1, 
 Real vx2, Real vy2, Real vz2) 
{
    Real r1_n = sqrt(sqr(x-x1) + sqr(y-y1) + sqr(z-z1));
    Real r2_n = sqrt(sqr(x-x2) + sqr(y-y2) + sqr(z-z2));

    Real v_r1 = (vx1 * (x - x1) + vy1 * (y - y1) + vz1 * (z - z1)) / r1_n;
    Real v_r2 = (vx2 * (x - x2) + vy2 * (y - y2) + vz2 * (z - z2)) / r2_n;

    return -(((dw12/w1 * speed_of_light + (dw12/w1-1) * v_r1) * r2_n + vy2 * (y - y2) + vz2 * (z - z2))/vx2 - x2);
}


template <typename Real = double>
auto F2(Real dw13, Real w1,
 Real x, Real y, Real z,
 Real x1,Real y1, Real z1, 
 Real x3,Real y3, Real z3, 
 Real vx1, Real vy1, Real vz1, 
 Real vx3, Real vy3, Real vz3) 
{
    Real r1_n = sqrt(sqr(x-x1) + sqr(y-y1) + sqr(z-z1));
    Real r3_n = sqrt(sqr(x-x3) + sqr(y-y3) + sqr(z-z3));

    Real v_r1 = (vx1 * (x - x1) + vy1 * (y - y1) + vz1 * (z - z1)) / r1_n;
    Real v_r3 = (vx3 * (x - x3) + vy3 * (y - y3) + vz3 * (z - z3)) / r3_n;

    return -(((dw13/w1 * speed_of_light + (dw13/w1-1) * v_r1) * r3_n + vx3 * (x - x3) + vz3 * (z - z3))/vy3 - y3);
}

template <typename Real = double>
auto F3(Real dw14, Real w1,
 Real x, Real y, Real z,
 Real x1,Real y1, Real z1, 
 Real x4,Real y4, Real z4, 
 Real vx1, Real vy1, Real vz1, 
 Real vx4, Real vy4, Real vz4) 
{
    Real r1_n = sqrt(sqr(x-x1) + sqr(y-y1) + sqr(z-z1));
    Real r4_n = sqrt(sqr(x-x4) + sqr(y-y4) + sqr(z-z4));

    Real v_r1 = (vx1 * (x - x1) + vy1 * (y - y1) + vz1 * (z - z1)) / r1_n;
    Real v_r4 = (vx4 * (x - x4) + vy4 * (y - y4) + vz4 * (z - z4)) / r4_n;

    return -(((dw14/w1 * speed_of_light + (dw14/w1-1) * v_r1) * r4_n + vy4 * (y - y4) + vx4 * (x - x4))/vz4 - z4);
}


template <typename Func1,typename Func2,typename Func3,typename Real>
auto zeidel_cycle3d(Real& x, Real& y, Real& z,Func1 f1, Func2 f2, Func3 f3)
{
    double x_ = x, y_ = y, z_ = z;

    x = f1(x_,y_,z_);
    y = f2(x_,y_,z_);
    z = f3(x_,y_,z_);
}

struct mat_dot
{
    boost::numeric::ublas::vector<double> pos, vel;

    mat_dot(double x = 0, double y = 0, double z = 0, double vx = 0,
            double vy = 0, double vz = 0)
            :pos(3), vel(3){
                pos(0) = x;
                pos(1) = y;
                pos(2) = z;
                vel(0) = vx;
                vel(1) = vy;
                vel(2) = vz;
            };
};

std::random_device rd;

auto generate_receiver(double R, double v)
{
    double angle = std::uniform_real_distribution<double>(0,2*3.14159265358979323)(rd);
    double vel_angle = std::uniform_real_distribution<double>(0,2*3.14159265358979323)(rd);

    double x = R*cos(angle);
    double y = R*sin(angle);
    double z = 0;
    double vx = v * sqrt(cos(vel_angle)*cos(vel_angle)/(1 + x*x/y/y));
    double vy = v * sqrt(cos(vel_angle)*cos(vel_angle)/(1 + y*y/x/x));
    double vz = v*sin(vel_angle);

    return mat_dot{x, y, z, vx, vy, vz};
}

auto generate_emitter(double R)
{
    double angle = std::uniform_real_distribution<double>(0,2*3.14159265358979323)(rd);
    return mat_dot{R*cos(angle),R*sin(angle), 0};
}

auto print(std::ostream& str, mat_dot& md)
{

    str << "pos:" << md.pos(0) << '\t' << md.pos(1) << '\t' << md.pos(2) << std::endl;
    str << "vel:" << md.vel(0) << '\t' << md.vel(1) << '\t' << md.vel(2) << std::endl;
}

int main()
{
    mat_dot receiver[4];
    receiver[0] = generate_receiver(1000,1);
    receiver[1] = generate_receiver(1000,1);
    receiver[2] = generate_receiver(700,1);
    receiver[3] = generate_receiver(700,1);

    mat_dot emitter = generate_emitter(100);

    print(std::cout,receiver[0]);
    print(std::cout,receiver[1]);
    print(std::cout,receiver[2]);
    print(std::cout,receiver[3]);

    std::cout << std::endl << "Emitter: " << std::endl;
    print(std::cout,emitter);
    std::cout << std::endl;

    double f = 1e9;


    double f1 = doppler_shift(f, receiver[0].vel, emitter.vel, receiver[0].pos, emitter.pos);
    double f2 = doppler_shift(f, receiver[1].vel, emitter.vel, receiver[1].pos, emitter.pos);
    double f3 = doppler_shift(f, receiver[2].vel, emitter.vel, receiver[2].pos, emitter.pos);
    double f4 = doppler_shift(f, receiver[3].vel, emitter.vel, receiver[3].pos, emitter.pos);

    double x{},y{},z{};

    for (auto i = 0; i < 100; ++i) {
        auto& r1 = receiver[0];
        auto& r2 = receiver[1];
        auto& r3 = receiver[2];
        auto& r4 = receiver[3];


        zeidel_cycle3d(x, y, z, 
        [&](double x, double y, double z)
        {return F1(f2-f1, f1, x, y, z, r1.pos(0), r1.pos(1),r1.pos(2), r2.pos(0), r2.pos(1),r2.pos(2), r1.vel(0), r1.vel(1),r1.vel(2), r2.vel(0), r2.vel(1),r2.vel(2));},

        [&](double x, double y, double z)
        {return F2(f3-f1, f1, x, y, z, r1.pos(0), r1.pos(1),r1.pos(2), r3.pos(0), r3.pos(1),r3.pos(2), r1.vel(0), r1.vel(1),r1.vel(2), r3.vel(0), r3.vel(1),r3.vel(2));},
        
        [&](double x, double y, double z)
        {return F3(f4-f1, f1, x, y, z, r1.pos(0), r1.pos(1),r1.pos(2), r4.pos(0), r4.pos(1),r4.pos(2), r1.vel(0), r1.vel(1),r1.vel(2), r4.vel(0), r4.vel(1),r4.vel(2));}
        );
    }
        std::cout << "result:" << x << ' ' << y << ' ' << z << std::endl;

}