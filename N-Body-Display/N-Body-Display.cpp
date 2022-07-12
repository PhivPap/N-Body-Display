#include <iostream>
#include "SFML/Graphics.hpp"
#include "util.h"
#include "body.h"
#include "quad.h"
#include "InputOutput\InputOutput.h"
#include "DisplayText.h"
#include <vector>
#include <chrono>
#include <omp.h>
#include <cmath>
#include <Windows.h>
#include <assert.h>

namespace KEY { 
    enum { ESC = 36, A = 0, D = 3, F = 5, G = 6, I = 8, O = 14, P = 15, R = 17, S = 18, W = 22, LESS = 49, MORE = 50, LEFT = 71, RIGHT = 72, UP = 73, DOWN = 74 };
}

const uint16_t MIN_BODY_RADIUS_PX = 1, MAX_BODY_RADIUS_PX = 255;
const double MIN_ITER_LEN = 1e-50, MAX_ITER_LEN = 1e50;

namespace FinalCFG { // never reset after init
    uint32_t h_res = 1600;
    uint32_t v_res = 900;
    double aspect_ratio = double(h_res) / v_res;
    const uint32_t max_fps = 60;
    const char* window_name = "N-Body";
    const char* font_path = "..\\data\\UbuntuMono-Regular.ttf";
    const uint32_t font_size = 26;
    const uint32_t font_y_spacing = 26;
    const double init_display_width = 1e14;
    const double max_display_width = 1e80;
    const double min_display_width = 1e-80;
    const double left_bound = -1e80;
    const double right_bound = 1e80;
    const double top_bound = -1e80;
    const double bottom_bound = 1e80;
    const double zoom_factor = 1.16;
    const double move_factor = 0.03;
    const uint32_t fps_update_freq = 30;     // fps info will be updated every 'fps_update_freq' frames

    std::string input_file = "../BodyFiles/in/a12a.tsv";
    std::string output_file = "../../../BodyFiles/out/bh_omp_out.tsv";
    uint32_t iterations = 100000;
    uint32_t thread_count = 8;
    double theta = 0.50;
    double theta_2 = theta * theta;
    uint32_t quad_pool_size = 0;

    void print(void) {
        std::cout << "\tWidth[px]: " << h_res << std::endl;
        std::cout << "\tHeight[px]: " << v_res << std::endl;
        std::cout << "\tAspect ratio: " << aspect_ratio << std::endl;
        std::cout << "\tMax FPS: " << max_fps << std::endl;
        std::cout << "\tWindow name: " << window_name << std::endl;
        std::cout << "\tInput file: " << input_file << std::endl;
        std::cout << "\tOutput file: " << output_file << std::endl;
        std::cout << "\tTotal iterations (steps): " << iterations << std::endl;
        std::cout << "\tRequested Threads (OpenMP): " << thread_count << std::endl;
        std::cout << "\tTheta: " << theta << std::endl;
    }
};

namespace MutCFG {
    bool camera_follow_body = false;
    uint32_t body_radius_px = 2;
    bool body_coloring = true;
    double iter_len = 1e4;

    void print(void) {
        std::cout << "\tBody radius[px]: " << body_radius_px << std::endl;
        const char* s = body_coloring ? "enabled" : "disabled";
        std::cout << "\tBody coloring: " << s << std::endl;
        std::cout << "\tTime step[s]: " << iter_len << std::endl;
    }

    namespace GridStatus {
        const uint8_t option_count = 3;
        enum GridStatus { None, XY, Standard };
        GridStatus grid_status = None;

        inline std::string grid_status_to_str(GridStatus gs) {
            switch (gs) {
            case None:
                return "None";
                break;
            case XY:
                return "X-Y axis";
                break;
            case Standard:
                return "Standard";
                break;
            default:
                assert(0);
                break;
            }
        }
    }
};



// ------------------------ SIMULATION ------------------------
const double G = 6.67e-11; // Gravitational constant
double merge_distance;
double max_abs_vel = 1e8;
uint32_t body_count_gui = 0;

void parse_input(const std::string& input_path, std::vector<Body>& bodies){
    try {
        IO::Parser parser(input_path);
        auto io_body = IO::Body();
        while (parser.next_body_info(io_body) == 0)
            bodies.push_back(Body(io_body.mass, { io_body.x, io_body.y }, io_body.vel_x, io_body.vel_y));
    }
    catch (const std::string& e){
        std::cout << "Parse error: " << e << std::endl;
        exit(1);
    }
}

void write_output(const std::string& output_path, const std::vector<Body>& bodies){
    try {
        IO::Writer writer(output_path);
        uint32_t id = 0;
        for (const auto& body : bodies)
            writer.write_body(IO::Body(id++, body.mass, body.coords.x, body.coords.y, body.vel_x, body.vel_y));
    }
    catch (const std::string& e){
        std::cout << "Write error: " << e << std::endl;
        exit(1);
    }
}

Area update_body_positions_get_area(Body* bodies, uint32_t body_count, double time_step){
    Area area(DBL_MAX, -DBL_MAX, DBL_MAX, -DBL_MAX);
    for (uint32_t i = 0; i < body_count; i++){
        Body& body = bodies[i];
        if (body.ignore)
            continue;
        auto delta_x = body.vel_x * time_step;
        auto delta_y = body.vel_y * time_step;
        body.coords.x += delta_x;
        body.coords.y += delta_y;
        area.x1 = std::min(area.x1, body.coords.x);
        area.x2 = std::max(area.x2, body.coords.x);
        area.y1 = std::min(area.y1, body.coords.y);
        area.y2 = std::max(area.y2, body.coords.y);
    }
    return area;
}

bool body_collision(Body* body1, Body* body2) {
    Body* big_body;
    Body* small_body;
    if (body1->mass > body2->mass){
        big_body = body1;
        small_body = body2;
    }
    else {
        big_body = body2;
        small_body = body1;
    }
    if (std::abs(small_body->vel_x - big_body->vel_x) < max_abs_vel && std::abs(small_body->vel_y - big_body->vel_y) < max_abs_vel)
        return false;
    double total_mass = big_body->mass + small_body->mass;
    big_body->mass = total_mass;
    big_body->vel_x = (small_body->mass * small_body->vel_x) / total_mass;
    big_body->vel_y = (small_body->mass * small_body->vel_y) / total_mass;
    small_body->ignore = true;
    body_count_gui--;
    return true;
}

void compute_body2body_attraction(Body* body1, Body* body2, double& Fx, double& Fy) {
    if (body1 == body2)
        return;
    const double distance_2 = body1->coords.distance_to_2(body2->coords);
    const double distance = std::sqrt(distance_2);
    /*if (distance < merge_distance) {
        if (!body1->ignore && !body2->ignore && body_collision(body1, body2))
            return;
    }*/
    const double F = (G * body1->mass * body2->mass) / distance_2;
    Fx += F * (body2->coords.x - body1->coords.x) / distance;
    Fy += F * (body2->coords.y - body1->coords.y) / distance;
}

void compute_body2quad_attraction(const Body* body, const Quad* quad, double& Fx, double& Fy, double distance_2) {
    const double distance = std::sqrt(distance_2);
    const double F = (G * body->mass * quad->mass) / distance_2;
    Fx += F * (quad->center_of_mass.x - body->coords.x) / distance;
    Fy += F * (quad->center_of_mass.y - body->coords.y) / distance;
}

void compute_body_forces(Quad* quad, Body* body, double& Fx, double& Fy) {
    const auto quad_body_count = quad->body_count;
    if (quad_body_count == 0)
        return;
    if (quad_body_count == 1){
        compute_body2body_attraction(body, quad->contained_bodies.front(), Fx, Fy);
        return;
    }

    // magic formula check https://en.wikipedia.org/wiki/Barnes%E2%80%93Hut_simulation
    const auto quad_diag_2 = quad->diag_len_2;
    const auto distance_2 = body->coords.distance_to_2(quad->center_of_mass);

    if (quad_diag_2 < FinalCFG::theta_2 * distance_2){
        compute_body2quad_attraction(body, quad, Fx, Fy, distance_2);
    }
    else {
        compute_body_forces(quad->top_left_quad, body, Fx, Fy);
        compute_body_forces(quad->top_right_quad, body, Fx, Fy);
        compute_body_forces(quad->bot_left_quad, body, Fx, Fy);
        compute_body_forces(quad->bot_right_quad, body, Fx, Fy);
    }
}

void update_body_velocities(Quad* root, Body* bodies, uint32_t body_count, double time_step) {
    #pragma omp parallel for default(none)              \
    firstprivate(root, bodies, body_count, time_step)   \
    schedule(dynamic)
    for (int32_t i = 0; i < body_count; i++){
        Body& body = bodies[i];
        if (body.ignore)
            continue;
        double Fx = 0, Fy = 0;
        compute_body_forces(root, &body, Fx, Fy);
        body.vel_x += (Fx / body.mass) * time_step;
        body.vel_y += (Fy / body.mass) * time_step;
    }
}

void simulate(Body* bodies, uint32_t body_count, double time_step){
    merge_distance = time_step * 1e6;   // the more time between time_steps, the more severe collision check becomes
    Quad::set_pool(FinalCFG::quad_pool_size);
    const auto area = update_body_positions_get_area(bodies, body_count, time_step);
    Quad root(bodies, body_count, area);
    update_body_velocities(&root, bodies, body_count, time_step);
    Quad::reset_pool();
}

// ------------------------- GRAPHICS -------------------------
uint32_t actual_threads = 1;
std::vector<DisplayText> display_texts;
double elapsed = 0.0;
double min_mass_log10, max_mass_log10;
double sim_elapsed = 0;
std::string iter_len_str;
std::string sim_elapsed_str;
std::string display_width_str, display_height_str;
std::string time_ratio_str;
std::string cpu_usage_str;
std::string grid_status_str;
uint32_t frames = 0;
double fps;
double universe_width = 0.0, universe_height = 0.0, universe_x1 = 0.0, universe_y1 = 0.0;
double universe_width_px = 0.0, universe_height_px = 0.0, universe_offset_x_px = 0.0, universe_offset_y_px = 0.0;
std::chrono::high_resolution_clock::time_point start;
sf::Font def_font;
sf::Rect<double> display_bound;

void init_gui(const Body* bodies, uint32_t body_count) {
    if (!def_font.loadFromFile(FinalCFG::font_path))
        exit(1);
    const uint32_t cfg_stats_dist = 200;
    float text_y = 0;
    body_count_gui = body_count;
    grid_status_str = MutCFG::GridStatus::grid_status_to_str(MutCFG::GridStatus::grid_status);

    display_texts.push_back(DisplayText("Iterations:  ", &FinalCFG::iterations, DisplayText::RefType::UINT32, { 5, text_y += FinalCFG::font_y_spacing }, def_font, FinalCFG::font_size));
    display_texts.push_back(DisplayText("Time step:   ", &iter_len_str, DisplayText::RefType::STR, { 5, text_y += FinalCFG::font_y_spacing }, def_font, FinalCFG::font_size));
    display_texts.push_back(DisplayText("Theta:       ", &FinalCFG::theta, DisplayText::RefType::DOUBLE, { 5, text_y += FinalCFG::font_y_spacing }, def_font, FinalCFG::font_size));
    display_texts.push_back(DisplayText("Threads:     ", &actual_threads, DisplayText::RefType::UINT32, { 5, text_y += FinalCFG::font_y_spacing }, def_font, FinalCFG::font_size));
    display_texts.push_back(DisplayText("Bodies:      ", &body_count_gui, DisplayText::RefType::UINT32, { 5, text_y += FinalCFG::font_y_spacing }, def_font, FinalCFG::font_size));
    display_texts.push_back(DisplayText("Max FPS:     ", &FinalCFG::max_fps, DisplayText::RefType::UINT32, { 5, text_y += FinalCFG::font_y_spacing }, def_font, FinalCFG::font_size));
    display_texts.push_back(DisplayText("Width[px]:   ", &FinalCFG::h_res, DisplayText::RefType::UINT32, { 5, text_y += FinalCFG::font_y_spacing }, def_font, FinalCFG::font_size));
    display_texts.push_back(DisplayText("Heigh[px]:   ", &FinalCFG::v_res, DisplayText::RefType::UINT32, { 5, text_y += FinalCFG::font_y_spacing }, def_font, FinalCFG::font_size));
    display_texts.push_back(DisplayText("Body R[px]:  ", &MutCFG::body_radius_px, DisplayText::RefType::UINT32, { 5, text_y += FinalCFG::font_y_spacing }, def_font, FinalCFG::font_size));
    display_texts.push_back(DisplayText("Colorcoding: ", &MutCFG::body_coloring, DisplayText::RefType::BOOL, { 5, text_y += FinalCFG::font_y_spacing }, def_font, FinalCFG::font_size));
    display_texts.push_back(DisplayText("Track heavy: ", &MutCFG::camera_follow_body, DisplayText::RefType::BOOL, { 5, text_y += FinalCFG::font_y_spacing }, def_font, FinalCFG::font_size));
    display_texts.push_back(DisplayText("Grid:        ", &grid_status_str, DisplayText::RefType::STR, { 5, text_y += FinalCFG::font_y_spacing }, def_font, FinalCFG::font_size));
    text_y += cfg_stats_dist;
    display_texts.push_back(DisplayText("Elapsed[s]:  ", &elapsed, DisplayText::RefType::DOUBLE, { 5, text_y += FinalCFG::font_y_spacing }, def_font, FinalCFG::font_size));
    display_texts.push_back(DisplayText("Simulated:   ", &sim_elapsed_str, DisplayText::RefType::STR, { 5, text_y += FinalCFG::font_y_spacing }, def_font, FinalCFG::font_size));
    display_texts.push_back(DisplayText("Frames:      ", &frames, DisplayText::RefType::UINT32, { 5, text_y += FinalCFG::font_y_spacing }, def_font, FinalCFG::font_size));
    display_texts.push_back(DisplayText("CPU load:    ", &cpu_usage_str, DisplayText::RefType::STR, { 5, text_y += FinalCFG::font_y_spacing }, def_font, FinalCFG::font_size));
    display_texts.push_back(DisplayText("FPS:         ", &fps, DisplayText::RefType::DOUBLE, { 5, text_y += FinalCFG::font_y_spacing }, def_font, FinalCFG::font_size));
    display_texts.push_back(DisplayText("Time ratio:  ", &time_ratio_str, DisplayText::RefType::STR, { 5, text_y += FinalCFG::font_y_spacing }, def_font, FinalCFG::font_size));
    display_texts.push_back(DisplayText("Width:       ", &display_width_str, DisplayText::RefType::STR, { 5, text_y += FinalCFG::font_y_spacing }, def_font, FinalCFG::font_size));
    display_texts.push_back(DisplayText("Height:      ", &display_height_str, DisplayText::RefType::STR, { 5, text_y += FinalCFG::font_y_spacing }, def_font, FinalCFG::font_size));

    const auto display_height = FinalCFG::init_display_width * ((double)FinalCFG::v_res / FinalCFG::h_res);
    display_bound = sf::Rect<double>(-FinalCFG::init_display_width / 2, -display_height / 2, FinalCFG::init_display_width, display_height);

    double min_mass = DBL_MAX, max_mass = -DBL_MAX;
    for (uint32_t i = 0; i < body_count; i++) {
        const double body_mass = bodies[i].mass;
        if (body_mass > max_mass)
            max_mass = body_mass;
        if (body_mass < max_mass)
            min_mass = body_mass;
    }
    min_mass_log10 = log10(min_mass);
    max_mass_log10 = log10(max_mass);

    #pragma omp parallel
    {
        #pragma omp single
        actual_threads = omp_get_num_threads();
    }
}

// from https://stackoverflow.com/questions/23143693/retrieving-cpu-load-percent-total-in-windows-with-c
static float calc_cpu_load(unsigned long long idleTicks, unsigned long long totalTicks) {
    static unsigned long long _previousTotalTicks = 0;
    static unsigned long long _previousIdleTicks = 0;
    unsigned long long totalTicksSinceLastTime = totalTicks - _previousTotalTicks;
    unsigned long long idleTicksSinceLastTime = idleTicks - _previousIdleTicks;
    float ret = 1.0f - ((totalTicksSinceLastTime > 0) ? ((float)idleTicksSinceLastTime) / totalTicksSinceLastTime : 0);
    _previousTotalTicks = totalTicks;
    _previousIdleTicks = idleTicks;
    return ret;
}

static unsigned long long file_time_to_int64(const FILETIME& ft){
    return (((unsigned long long)(ft.dwHighDateTime)) << 32) | ((unsigned long long)ft.dwLowDateTime);
}

static double log(double x, double base) {
    return log(x) / log(base);
}

float get_cpu_load(void){
    FILETIME idleTime, kernelTime, userTime;
    return GetSystemTimes(&idleTime, &kernelTime, &userTime) ? calc_cpu_load(file_time_to_int64(idleTime), file_time_to_int64(kernelTime) + file_time_to_int64(userTime)) : -1.0f;
}

std::string seconds_to_string(double s) {
    if (s < 1e-6)       return DisplayText::double_to_str(s * 1e9) + " ns";
    if (s < 1e-3)       return DisplayText::double_to_str(s * 1e6) + " us";
    if (s < 1.0)        return DisplayText::double_to_str(s * 1e3) + " ms";
    if (s < 60.0)       return DisplayText::double_to_str(s) + " s";
    if (s < 3600.0)     return DisplayText::double_to_str(s / 60) + " min";
    if (s < 86400.0)    return DisplayText::double_to_str(s / 3600) + " hours";
    if (s < 31556952.0) return DisplayText::double_to_str(s / 86400) + " days";
    else                return DisplayText::double_to_str(s / 31556952) + " years";
}

std::string meters_to_string(double m) {
    if (m < 1e-6)               return DisplayText::double_to_str(m * 1e9) + " nm";
    if (m < 1e-3)               return DisplayText::double_to_str(m * 1e6) + " um";
    if (m < 1)                  return DisplayText::double_to_str(m * 1e3) + " mm";
    if (m < 1e3)                return DisplayText::double_to_str(m) + " m";
    if (m < 1495978707e2)       return DisplayText::double_to_str(m * 1e-3) + " km";
    if (m < 94607304725808e2)   return DisplayText::double_to_str(m / 1495978707e2) + " AU";
    else                        return DisplayText::double_to_str(m / 94607304725808e2) + " LY";
}

sf::Color mass_to_color(double kg) {
    static const double mass_log10_amplitude = max_mass_log10 - min_mass_log10;
    const double importance = (log10(kg) - min_mass_log10) / mass_log10_amplitude;
    const double importance_sup = 1 - importance;
    uint8_t red = 255;
    uint8_t green = 255 - importance_sup * 125;
    uint8_t blue = 255 - importance_sup * 125;
    return sf::Color(red, green, blue);
}

void display_rebound(void) {
    const auto tmp1 = FinalCFG::right_bound - display_bound.width;
    const auto tmp2 = FinalCFG::bottom_bound - display_bound.height;
    if (display_bound.left < FinalCFG::left_bound)
        display_bound.left = FinalCFG::left_bound;
    else if (display_bound.left > tmp1)
        display_bound.left = tmp1;
    if (display_bound.top < FinalCFG::top_bound)
        display_bound.top = FinalCFG::top_bound;
    else if (display_bound.top > tmp2)
        display_bound.top = tmp2;
}

void display_follow_body(const Body* follow) {
    display_bound.left = follow->coords.x - (display_bound.width / 2);
    display_bound.top = follow->coords.y - (display_bound.height / 2);
}

void move_display(int32_t dx, int32_t dy) {
    display_bound.left += dx * FinalCFG::move_factor * display_bound.width;
    display_bound.top += dy * FinalCFG::move_factor * display_bound.height;
}

void reset_display(void) {
    display_bound.left = -display_bound.width / 2;
    display_bound.top = -display_bound.height / 2;
}

void zoom_display(int32_t zooms) {
    if (zooms > 0 && display_bound.width >= FinalCFG::max_display_width)
        return;
    else if (display_bound.width <= FinalCFG::min_display_width)
        return;
    const auto resize = pow(FinalCFG::zoom_factor, zooms);
    const auto old_width = display_bound.width;
    const auto old_height = display_bound.height;
    display_bound.width *= resize;
    display_bound.height = display_bound.width * (1 / FinalCFG::aspect_ratio);
    display_bound.left += (old_width - display_bound.width) / 2;
    display_bound.top += (old_height - display_bound.height) / 2;
}

void update_counters(uint32_t iteration) {
    const auto now = std::chrono::high_resolution_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(now - start).count() / 1e9;

    iter_len_str = seconds_to_string(MutCFG::iter_len);

    sim_elapsed += MutCFG::iter_len;
    sim_elapsed_str = seconds_to_string(sim_elapsed);

    display_width_str = meters_to_string(display_bound.width);
    display_height_str = meters_to_string(display_bound.height);

    if (iteration % FinalCFG::fps_update_freq == 0) {
        static auto x_frames_ago = now;
        const auto x_frames_elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(now - x_frames_ago).count() / 1e9;
        fps = FinalCFG::fps_update_freq / x_frames_elapsed;
        time_ratio_str = seconds_to_string(MutCFG::iter_len * fps) + "/s";
        x_frames_ago = now;
        cpu_usage_str = DisplayText::double_to_str(get_cpu_load() * 100.0) + "%";
    }
}

void draw_standard_grid(sf::RenderWindow& window) {
    auto horizontal_line = sf::RectangleShape({ float(FinalCFG::h_res), 1 });
    auto vertical_line = sf::RectangleShape({ 1, float(FinalCFG::v_res) });
    horizontal_line.setFillColor(sf::Color(50, 50, 50));
    vertical_line.setFillColor(sf::Color(50, 50, 50));

    const double max_dim = std::min(display_bound.width, display_bound.height);
    const double spacing = pow(10, floor(log10(max_dim)) - 1);
    const double x2 = display_bound.width + display_bound.left;
    double x = ceil(display_bound.left / spacing) * spacing;
    while (x < x2) {
        const float x_ratio = (x - display_bound.left) / display_bound.width;
        vertical_line.setPosition({ x_ratio * FinalCFG::h_res, 0 });
        window.draw(vertical_line);
        x += spacing;
    }

    const double y2 = display_bound.height + display_bound.top;
    double y = ceil(display_bound.top / spacing) * spacing;
    while (y < y2) {
        const float y_ratio = (y - display_bound.top) / display_bound.height;
        horizontal_line.setPosition({ 0, y_ratio * FinalCFG::v_res });
        window.draw(horizontal_line);
        y += spacing;
    }
}

//void draw_standard_grid2(sf::RenderWindow& window) {
//    static const double spacing_factor = 5.0;
//    auto horizontal_line = sf::RectangleShape({ float(FinalCFG::h_res), 1 });
//    auto vertical_line = sf::RectangleShape({ 1, float(FinalCFG::v_res) });
//    horizontal_line.setFillColor(sf::Color(50, 50, 50));
//    vertical_line.setFillColor(sf::Color(50, 50, 50));
//    const double max_dim = std::min(display_bound.width, display_bound.height);
//    const double spacing = pow(spacing_factor, floor(log(max_dim, spacing_factor)) - 1);
//    const float pixel_spacing = spacing / max_dim;
//    double x = ceil(display_bound.left / spacing) * spacing;
//    double y = ceil(display_bound.top / spacing) * spacing;
//    const double x2 = display_bound.width + display_bound.left;
//    const double y2 = display_bound.height + display_bound.top;
//    float x_ratio = (x - display_bound.left) / display_bound.width;
//
//
//    while (x < x2) {
//        //const float x_ratio = (x - display_bound.left) / display_bound.width;
//        vertical_line.setPosition({ x_ratio * FinalCFG::h_res, 0 });
//        window.draw(vertical_line);
//        x_ratio += pixel_spacing;
//        x += spacing;
//    }
//
//
//
//}



void draw_xy_axis(sf::RenderWindow& window) {
    auto axis_x = sf::RectangleShape({ float(FinalCFG::h_res), 1 });
    auto axis_y = sf::RectangleShape({ 1, float(FinalCFG::v_res) });
    auto y = 0;
    auto x = 0;

    const float x_ratio = (x - display_bound.left) / display_bound.width;
    const float y_ratio = (y - display_bound.top) / display_bound.height;
    axis_x.setPosition({ 0, y_ratio * FinalCFG::v_res });
    axis_y.setPosition({ x_ratio * FinalCFG::h_res, 0 });
    axis_x.setFillColor(sf::Color(50, 50, 50));
    axis_y.setFillColor(sf::Color(50, 50, 50));
    window.draw(axis_x);
    window.draw(axis_y);
}

void draw_grid(sf::RenderWindow& window) {
    using namespace MutCFG::GridStatus;
    switch (grid_status) {
    case XY:
        draw_xy_axis(window);
        break;
    case Standard:
        draw_standard_grid(window);
        break;
    default:
        break;
    }
}

void draw_texts(sf::RenderWindow& window) {
    for (auto& dt : display_texts)
        dt.draw_to(window);
}

/* returns most massive object printed */
const Body* draw_bodies2(sf::RenderWindow& window, Body* bodies, uint32_t body_count) {
    double max_mass = -DBL_MAX;
    const Body* heaviest_body = nullptr;
    sf::CircleShape circle;
    circle.setRadius(MutCFG::body_radius_px);

    for (uint32_t i = 0; i < body_count; i++) {
        const Body& body = bodies[i];
        if (body.ignore)
            continue;
        if (!display_bound.contains(body.coords.x, body.coords.y))
            continue;
        if (body.mass > max_mass) {
            max_mass = body.mass;
            heaviest_body = &body;
        }
        const float x_ratio = (body.coords.x - display_bound.left) / display_bound.width;
        const float y_ratio = (body.coords.y - display_bound.top) / display_bound.height;
        circle.setPosition(x_ratio * FinalCFG::h_res, y_ratio * FinalCFG::v_res);
        const auto body_color = MutCFG::body_coloring ? mass_to_color(body.mass) : sf::Color::White;
        circle.setFillColor(body_color);
        window.draw(circle);
    }
    return heaviest_body;
}


void engine_loop(sf::RenderWindow& window, Body* bodies, uint32_t body_count) {
    sf::Event event;
    bool paused = false;
    const Body* follow = nullptr;
    for (uint32_t it = 0; it < FinalCFG::iterations; it++){
        if (!window.isOpen()) {
            std::cerr << "Internal error: window" << std::endl;
            exit(1);
        }
        int32_t zoom = 0;
        int32_t dx = 0, dy = 0;
        while (window.pollEvent(event)){
            switch (event.type) {
            case sf::Event::Closed:
                window.close();
                return;
            case sf::Event::MouseWheelScrolled:
                zoom = -event.mouseWheelScroll.delta;
                break;
            case sf::Event::KeyPressed:
                //std::cout << event.key.code << std::endl;
                switch (event.key.code){
                case KEY::ESC:
                    if (paused) {
                        paused = false;
                    }
                    else {
                        DisplayText("Paused", nullptr, DisplayText::RefType::NONE, { (float(FinalCFG::h_res / 2)) - 60 , 50 }, def_font, 40).draw_to(window);
                        window.display();
                        paused = true;
                    }
                    break;
                case KEY::F:
                    MutCFG::camera_follow_body = !MutCFG::camera_follow_body;
                    break;
                case KEY::G: {
                    using namespace MutCFG::GridStatus;
                    grid_status = static_cast<GridStatus>((grid_status + 1) % option_count);
                    grid_status_str = grid_status_to_str(grid_status);
                    break;
                }
                case KEY::O:
                    if (MutCFG::body_radius_px > MIN_BODY_RADIUS_PX)
                        MutCFG::body_radius_px--;
                    break;
                case KEY::P:
                    if (MutCFG::body_radius_px < MAX_BODY_RADIUS_PX)
                        MutCFG::body_radius_px++;
                    break;
                case KEY::I:
                    MutCFG::body_coloring = !MutCFG::body_coloring;
                    break;
                case KEY::LESS:
                    if (MutCFG::iter_len > MIN_ITER_LEN)
                        MutCFG::iter_len /= 2;
                    break;
                case KEY::MORE:
                    if (MutCFG::iter_len < MAX_ITER_LEN)
                        MutCFG::iter_len *= 2;
                    break;
                case KEY::UP:
                case KEY::W:
                    dy--;
                    break;
                case KEY::DOWN:
                case KEY::S:
                    dy++;
                    break;
                case KEY::LEFT:
                case KEY::A:
                    dx--;
                     break;
                case KEY::RIGHT:
                case KEY::D:
                    dx++;
                    break;
                case KEY::R:
                    dx = dy = 0;
                    if (!MutCFG::camera_follow_body)
                        reset_display();
                }
                break;
            }
        }

        if (paused) {
            it--;
            continue;
        }
        window.clear();
        if (zoom != 0)
            zoom_display(zoom);
        if (MutCFG::camera_follow_body && follow != nullptr)
            display_follow_body(follow);
        else if (dx != 0 || dy != 0)
            move_display(dx, dy);
        display_rebound();
        update_counters(it);
        follow = draw_bodies2(window, bodies, body_count);
        draw_grid(window);
        draw_texts(window);
        window.display();
        simulate(bodies, body_count, MutCFG::iter_len);
        frames++;
    }
    window.clear();
    draw_bodies2(window, bodies, body_count);
    draw_texts(window);
    DisplayText("Simulation Complete", nullptr, DisplayText::RefType::NONE, { (float(FinalCFG::h_res / 2)) - 190 , 50 }, def_font, 40).draw_to(window);
    window.display();


    while (window.isOpen()) {
        while (window.pollEvent(event)){
            switch (event.type) {
            case sf::Event::Closed:
                window.close();
                return;
            }
        }
    }
}

void parse_args(int argc, const char** argv){
    int32_t idx;
    IO::ArgParser arg_parser(argc, argv);
    try{
        if ((idx = arg_parser.get_next_idx("-in")) > 0)
            FinalCFG::input_file = arg_parser.get(idx);

        if ((idx = arg_parser.get_next_idx("-out")) > 0)
            FinalCFG::output_file = arg_parser.get(idx);

        if ((idx = arg_parser.get_next_idx("-it")) > 0)
            FinalCFG::iterations = std::stoi(arg_parser.get(idx));

        if ((idx = arg_parser.get_next_idx("-it_len")) > 0)
            MutCFG::iter_len = std::stod(arg_parser.get(idx));

        if ((idx = arg_parser.get_next_idx("-theta")) > 0){
            FinalCFG::theta = std::stod(arg_parser.get(idx));
            FinalCFG::theta_2 = FinalCFG::theta * FinalCFG::theta;
        }

        if ((idx = arg_parser.get_next_idx("-threads")) > 0)
            FinalCFG::thread_count = std::stod(arg_parser.get(idx));

        if ((idx = arg_parser.get_next_idx("-h_res")) > 0) {
            using namespace FinalCFG;
            h_res = std::stoi(arg_parser.get(idx));
            aspect_ratio = double(h_res) / v_res;
        }

        if ((idx = arg_parser.get_next_idx("-v_res")) > 0) {
            using namespace FinalCFG;
            v_res = std::stoi(arg_parser.get(idx));
            aspect_ratio = double(h_res) / v_res;
        }

        if ((idx = arg_parser.get_next_idx("-max_fps")) > 0)
            FinalCFG::h_res = std::stoi(arg_parser.get(idx));

    }
    catch (const std::string& ex){
        std::cout << "Argument parsing exception: " << ex << std::endl;
        exit(1);
    }
    catch (const std::invalid_argument& e){
        std::cout << "Argument parsing exception: invalid argument." << std::endl;
        exit(1);
    }
    catch (const std::out_of_range& e){
        std::cout << "Argument parsing exception: out of range." << std::endl;
        exit(1);
    }
}

void print_info(void) {
    std::cout << "Configuration:\n";
    FinalCFG::print();
    MutCFG::print();    
    std::cout << "\nKeybinds:\n";
    std::cout << "\tO: decrease body radius [px]\n";
    std::cout << "\tP: increase body radius [px]\n";
    std::cout << "\tI: enable or disable body coloring\n";
    std::cout << "\t<scroll up>: zoom in\n";
    std::cout << "\t<scroll down>: zoom out\n";
    std::cout << "\t<: slowdown the simulation (halves timestep)\n";
    std::cout << "\t>: speedup the simulation (doubles timestep)\n";
    std::cout << "\tW/<up arrow>: moves view by a small factor north\n";
    std::cout << "\tS/<down arrow>: moves view by a small factor south\n";
    std::cout << "\tA/<left arrow>: moves view by a small factor west\n";
    std::cout << "\tD/<right arrow>: moves view by a small factor east\n";
    std::cout << "\tR: Resets view to point at (x, y) = (0, 0)\n";
    std::cout << "\tF: Camera follows heaviest body in sight\n";
    std::cout << "\tG: Toggle grid\n";

}

int main(int argc, const char** argv){
    std::vector<Body> bodies;
    parse_args(argc, argv);
    print_info();
    omp_set_num_threads(FinalCFG::thread_count);
    parse_input(FinalCFG::input_file, bodies);
    FinalCFG::quad_pool_size = 400 + bodies.size() * 4;
    init_gui(bodies.data(), bodies.size());
    sf::RenderWindow window(sf::VideoMode(FinalCFG::h_res, FinalCFG::v_res), FinalCFG::window_name);
    window.setFramerateLimit(FinalCFG::max_fps);
    start = std::chrono::high_resolution_clock::now();
    engine_loop(window, bodies.data(), bodies.size());
    const auto end = std::chrono::high_resolution_clock::now();
    const auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1e9;
    //std::cout << "Seconds: " << elapsed << std::endl;
    //std::cout << "Seconds per iteration: " << elapsed / FinalCFG::iterations << std::endl;
    write_output(FinalCFG::output_file, bodies);
    std::cout << "Jobs done!" << std::endl;
}