#include <Common.hpp>
#include <Point2D.hpp>
#include <cmath>
#include <cstdlib>
#include <iostream>

#include <Vizior/Vizior.hpp>
#include <Vizior/Window.hpp>
#include <Vizior/ImageBuilder.hpp>
#include <memory>

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

std::shared_ptr<vzr::Window> win;
std::shared_ptr<vzr::ImageBuilder> imb;


double start_angle = 3.141592/16;
double gravity = 9.81;

vzr::Point2D center{800, 400};
double radius = 250;
vzr::Point2D absolute_pos{800.0, 400.0+radius};
double speed = 100;
vzr::Vec2 velocity = (absolute_pos - center).rotatedBy(-vzr::C_PI/2).normalized()*speed;
double attraction_strength = speed*speed / radius;

float timestep = 1.0/60.0;

int steps_per_frame = 1;


void update_forward_euler(){
  //
  // xn+1 = xn + dt*vn
  //
  // vn+1 = vn + dt*fn/m
  //

  for(int i = 0; i < steps_per_frame; i++){
    vzr::Vec2 strength = (center - absolute_pos).normalized() * attraction_strength;
  
    vzr::Vec2 next_vel = velocity + strength*double(timestep);

    vzr::Vec2 pos = absolute_pos - center;
  
    vzr::Vec2 next_pos = pos + velocity*double(timestep); // Need to add ALL the operators needed to avoid this...
  
  
    absolute_pos = center + next_pos;
    velocity = next_vel;  
  }
}

void update_symplectic_euler(){
  //
  // xn+1 = xn + dt*vn+1
  //
  // vn+1 = vn + dt*fn/m
  //

  for(int i = 0; i < steps_per_frame; i++){
    vzr::Vec2 strength = (center - absolute_pos).normalized() * attraction_strength;
  
    vzr::Vec2 next_vel = velocity + timestep*strength;

    vzr::Vec2 pos = absolute_pos - center;
  
    vzr::Vec2 next_pos = pos + timestep*next_vel;
  
    absolute_pos = center + next_pos;
    velocity = next_vel;  
  }
}

void inverse_22_mat(double *a, double *b, double *c, double *d){
  double det = (*a)*(*d) - (*b)*(*c);
  double new_a = *d/det;
  double new_b = -*b/det;
  double new_c = -*c/det;
  double new_d = *a/det;

  *a = new_a;
  *b = new_b;
  *c = new_c;
  *d = new_d;
}

void update_backwards_euler(){
  //
  // xn+1 = xn + dt*vn+1
  //
  // vn+1 = vn + dt*fn+1/m
  //
  // Since xn+1 depends on vn+1, which depends on fn+1 which is a function of x, we need to solve a system/non-linear equation
  //
  // We can work the equaton to be 0 = xn - xn+1 + vn*dt + dt^2 * Fn+1 = g(xn+1)
  // Remember that g: R^2 -> R^2 and xn+1 = (x,y)
  // Also Fn+1 = -k * (x/||X||, y/||X||) (capital X being the vector and lowercase x being the coordinate in this specific line)
  // 
  // We need the jabobian, which we can easily compute
  // dg1/dx = 1 + y^2/((x^2+y^2)^(3/2))
  // dg2/dx =     -xy/((x^2+y^2)^(3/2))
  // dg1/dy =     -xy/((x^2+y^2)^(3/2))
  // dg2/dy = 1 + x^2/((x^2+y^2)^(3/2))
  //
  // Note that in our (very simple) case, we take the inverse of the matrix
  // In usual cases we would solve the system J(xn+1 - xn) = F using other methods since computing the inverse vecomes
  // very expensive for very large matrices

  for(int i = 0; i < steps_per_frame; i++){
    vzr::Vec2 pos = absolute_pos - center;
    
    // Newton's method
    vzr::Vec2 guess = pos;
    vzr::Vec2 strength = (-guess).normalized() * attraction_strength;
    double dt2 = timestep*timestep;
    auto g = [&pos, &strength,  dt2](vzr::Vec2& xi)-> vzr::Vec2{
      return pos - xi + velocity*timestep + strength*dt2;
    };
    vzr::Vec2 gGuess = g(guess);
    int iter = 0;
    while( gGuess.length() > 0.0000001){
      std::cout << "Newton from new guess " << guess.x << "," << guess.y << " gGuess is " << gGuess.x << " " << gGuess.y << std::endl;
      double det = guess.length() * guess.lengthSquared();
      double a, dg1dx = - 1 - dt2*attraction_strength*guess.y*guess.y/det; a = dg1dx;
      double b, dg1dy =     + dt2*attraction_strength*guess.y*guess.x/det; b = dg1dy;
      double c, dg2dx =     + dt2*attraction_strength*guess.y*guess.x/det; c = dg2dx;
      double d, dg2dy = - 1 - dt2*attraction_strength*guess.x*guess.x/det; d = dg2dy;

      inverse_22_mat(&a, &b, &c, &d);

      // Manual mat mult
      guess.x -= gGuess.x * a + gGuess.y * b;
      guess.y -= gGuess.x * c + gGuess.y * d;
      
      strength = (-guess).normalized() * attraction_strength;
      gGuess = g(guess);
      iter++;
    }
    std::cout << "Finished newton in " << iter << " iterations guess is " << guess.x << "," << guess.y << std::endl; 

    vzr::Vec2 next_pos = guess;

  
    vzr::Vec2 next_vel = (next_pos - pos) / timestep;

  
  
    absolute_pos = center + next_pos;
    velocity = next_vel;  
  }
}

void main_loop(int frame_count, double _d){
  ImGui_ImplOpenGL3_NewFrame();
  ImGui_ImplGlfw_NewFrame();
  ImGui::NewFrame();

  ImGui::Begin("Simulation parameters");
  ImGui::Text("Frame number %d", frame_count);
  ImGui::SliderFloat("Timestep", &timestep, 1.0/600.0, 1.0);
  ImGui::SliderInt("Steps per frame", &steps_per_frame, 0, 100);
  ImGui::End();

  update_backwards_euler();

  imb->drawCircleOutline(center, radius, 2, vzr::BLK);
  imb->drawCircle(absolute_pos, 20, {255,0,0,128});
  imb->drawCircle(absolute_pos, 2, {255,255,255,128});
  imb->drawCircle(center, 5, vzr::BLK);
}

int main(){

  
  
  win = std::make_shared<vzr::Window>(1600,800,"LinkTest");
  imb = std::make_shared<vzr::ImageBuilder>();
  imb->setBackgroundColor({210, 180, 180, 255});
  win->setSource(imb);
  vzr::registerWindow(win);
  vzr::setLoopFunc(main_loop);
  
  IMGUI_CHECKVERSION();
  ImGui::CreateContext();

  // Setup Platform/Renderer backends
  ImGui_ImplGlfw_InitForOpenGL(win->getGLFWWindow(), true);          // Second param install_callback=true will install GLFW callbacks and chain to existing ones.
  ImGui_ImplOpenGL3_Init();
  
  win->set_pre_swapbuffer_callback([]{
    ImGui::Render();
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
  });

  vzr::Start();

  ImGui_ImplOpenGL3_Shutdown();
  ImGui_ImplGlfw_Shutdown();
  ImGui::DestroyContext();
}
