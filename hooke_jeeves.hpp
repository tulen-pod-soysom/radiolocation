
#include <algorithm>
#include <iterator>


template <typename Func, typename Vec>
auto hooke_jeeves_optimizer(Func f, Vec starting_point, unsigned dimensions,double precision = 1e-9,double starting_step = 1, double speed_factor = 2.0)
{
  Vec steps(dimensions,starting_step);
  auto& point1 = starting_point;
  auto  point2 = starting_point;
  auto  point3 = starting_point;

  do
  {
    // Исследующий поиск
    for (auto i = 0; i < dimensions; ++i) {
      Vec point_plus = point2;
      Vec point_minus = point2;

      *(std::begin(point_plus) + i) += *(std::begin(steps) + i);
      *(std::begin(point_minus) + i) -= *(std::begin(steps) + i);

      auto f_plus = f(point_plus);
      auto f_middle = f(point2);
      auto f_minus = f(point_minus);

      if ((f_middle <= f_plus) && (f_middle <= f_minus)) {*(std::begin(steps) + i) /= 2.0;}
      else if (f_plus < f_minus) {point2 = point_plus;}
      else {point2 = point_minus;}
    }

    // Поиск по образцу (шаблону)
    point3 = point1;

    for (auto i = 0; i < dimensions; ++i)
    {
      auto& c1 = *(std::begin(point1) + i);
      auto& c2 = *(std::begin(point2) + i);
      auto& c3 = *(std::begin(point3) + i);
      c3 += speed_factor*(c2 - c1);
    }

    for (auto i = 0; i < dimensions; ++i) {
      Vec point_plus = point3;
      Vec point_minus = point3;

      *(std::begin(point_plus) + i) += *(std::begin(steps) + i);
      *(std::begin(point_minus) + i) -= *(std::begin(steps) + i);

      auto f_plus = f(point_plus);
      auto f_middle = f(point3);
      auto f_minus = f(point_minus);

      if ((f_middle <= f_plus) && (f_middle <= f_minus)) {point1 = point2;}
      else if (f_plus < f_minus) {point1 = point2; point2 = point_plus;}
      else {point1 = point2; point2 = point_minus;}
    }
  }
  while (std::any_of(steps.begin(),steps.end(),[&](auto& a){return a > precision;}));

  return point2;
}

template <typename Func, typename Vec>
auto hooke_jeeves_optimizer_simple(Func f, Vec starting_point, unsigned dimensions,double precision = 1e-9,double starting_step = 1)
{
  Vec steps(dimensions,starting_step);
  auto& point1 = starting_point;


  do
    for (auto i = 0; i < dimensions; ++i) {
      Vec point_plus = point1;
      Vec point_minus = point1;

      *(std::begin(point_plus) + i) += *(std::begin(steps) + i);
      *(std::begin(point_minus) + i) -= *(std::begin(steps) + i);

      auto f_plus = f(point_plus);
      auto f_middle = f(point1);
      auto f_minus = f(point_minus);

      if ((f_middle <= f_plus) && (f_middle <= f_minus)) {*(std::begin(steps) + i) /= 2.0;}
      else if (f_plus < f_minus) {point1 = point_plus;}
      else {point1 = point_minus;}
    }
    while (std::any_of(steps.begin(),steps.end(),[&](auto& a){return a > precision;}));

    return point1;
}
