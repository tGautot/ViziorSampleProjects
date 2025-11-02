#include <cmath>
#include <iostream>

#include <Vizior/Vizior.hpp>
#include <Vizior/Window.hpp>
#include <Vizior/ImageBuilder.hpp>
#include <memory>

std::shared_ptr<vzr::Window> win;
std::shared_ptr<vzr::ImageBuilder> imb;


double start_angle = 3.141592/16;
double gravity = 9.81;

vzr::Point2D attach_point1{400, 600};
vzr::Point2D attach_point2{1200, 600};
double attach_len = 400;

double total_time = 0;


void loop_formula(int frame_count, double timestep){
  total_time += timestep;
  double current_angle = start_angle * std::cos(std::sqrt(gravity/4) * total_time);
  vzr::Point2D ball_pos = attach_point1;
  ball_pos.y -= attach_len;
  ball_pos.rotateAroundBy(attach_point1, current_angle);

  vzr::Color col = {20,30,30,255};
  imb->drawLine({attach_point1, ball_pos}, 5, col, false);
  imb->drawCircle(ball_pos, 50, col);
}

double current_radial_velocity = 0.0;

void loop_simulation(int frame_count, double timestep){
  static double current_angle = start_angle;
  double radial_accel = -gravity * std::sin(current_angle) / 4;
  current_radial_velocity += radial_accel*timestep;
  current_angle += current_radial_velocity*timestep;

  vzr::Point2D ball_pos = attach_point2;
  ball_pos.y -= attach_len;
  ball_pos.rotateAroundBy(attach_point2, current_angle);

  vzr::Color col = {20,30,30,255};
  //std::cout << "Angle of second ball is " << current_angle << " final pos is {" << ball_pos.x << "," << ball_pos.y <<"}" << std::endl;
  imb->drawLine({attach_point2, ball_pos}, 5, col, false);
  imb->drawCircle(ball_pos, 50, col);
}

void main_loop(int frame_count, double timestep){
  loop_formula(frame_count, timestep);
  loop_simulation(frame_count, timestep);
}

int main(){
  win = std::make_shared<vzr::Window>(1600,800,"LinkTest");
  imb = std::make_shared<vzr::ImageBuilder>();
  imb->setBackgroundColor({210, 180, 180, 255});
  win->setSource(imb);
  vzr::registerWindow(win);
  vzr::setLoopFunc(main_loop);
  vzr::Start();
}
