#include "quad.h"
#include <assert.h>
#include <chrono>

const uint8_t MAX_DEPTH = 40;
//std::vector<Quad> Quad::pool;
Quad *Quad::pool;
uint32_t Quad::pool_size;
uint32_t Quad::pool_idx;

Quad::Quad(void) : mass(0), body_count(0), center_of_mass({0, 0}){
    contained_bodies.reserve(40);
}

// this constructor is used to generate the root of the quad tree
Quad::Quad(Body *bodies, uint32_t body_count, const Area &area) : center_of_mass({0, 0}), body_count(0),
                                                                  mass(0), area(area) {
    diag_len_2 = area.diagonal_length_2();
    contained_bodies.reserve(body_count);
    for (uint32_t i = 0; i < body_count; i++) {
        Body* b = &(bodies[i]);
        if (!b->ignore)
            insert_body(b);
    }
    compute_bhtree_recursive(0);
}

Quad::~Quad() {}

void Quad::set_area(const Area &area){
    this->area = area;
    diag_len_2 = area.diagonal_length_2();
}

void Quad::insert_body(Body *body){
    contained_bodies.push_back(body);
    mass += body->mass;
    body_count++;
}

void Quad::compute_bhtree_recursive(uint8_t depth){
    if (body_count == 0)
        return;
    else if (body_count == 1){
        center_of_mass = contained_bodies[0]->coords;
        return;
    }

    if (depth > MAX_DEPTH) {
        body_count = 1;
        printf("Tree Depth reached %d, dropping %llu bodies from next force calc.\n", MAX_DEPTH, contained_bodies.size());
        return;
    }

    const Point &center = area.get_center();

    auto my_pool_index = Quad::pool_get_idx();
    auto subtree = &(Quad::pool[my_pool_index]);

    //top_left_idx = my_pool_index;

    top_left_quad = subtree;
    top_right_quad = &(subtree[1]);
    bot_left_quad = &(subtree[2]);
    bot_right_quad = &(subtree[3]);

    top_left_quad->set_area({area.x1, center.x, area.y1, center.y});
    top_right_quad->set_area({center.x, area.x2, area.y1, center.y});
    bot_left_quad->set_area({area.x1, center.x, center.y, area.y2});
    bot_right_quad->set_area({center.x, area.x2, center.y, area.y2});

    const uint32_t n = contained_bodies.size();
    auto cb = contained_bodies.data();

    for (uint32_t i = 0; i < n; i++) {
        Body *body = cb[i];
        const Point &coords = body->coords;

        if (coords.x > center.x){
            if (coords.y > center.y)
                bot_right_quad->insert_body(body);
            else
                top_right_quad->insert_body(body);
        }
        else {
            if (coords.y > center.y)
                bot_left_quad->insert_body(body);
            else
                top_left_quad->insert_body(body);
        }
    }

    depth++;
    top_left_quad->compute_bhtree_recursive(depth);
    top_right_quad->compute_bhtree_recursive(depth);
    bot_left_quad->compute_bhtree_recursive(depth);
    bot_right_quad->compute_bhtree_recursive(depth);

    center_of_mass = Point::get_center_of_mass(
        Point::get_center_of_mass(
            top_left_quad->center_of_mass, top_left_quad->mass,
            top_right_quad->center_of_mass, top_right_quad->mass),
        top_left_quad->mass + top_right_quad->mass,
        Point::get_center_of_mass(
            bot_left_quad->center_of_mass, bot_left_quad->mass,
            bot_right_quad->center_of_mass, bot_right_quad->mass),
        bot_left_quad->mass + bot_right_quad->mass);
}



void Quad::set_pool(uint32_t init_size){
    /*pool.reserve(init_size);*/
    pool_idx = 0;
    pool_size = init_size;
    pool = new Quad[pool_size];
}

void Quad::reset_pool(void){
    delete[] pool;
}

uint32_t Quad::pool_get_idx(void){
    const uint32_t curr_idx = pool_idx;
    pool_idx += 4;
    assert(pool_idx < pool_size);
    return curr_idx;
    /*uint32_t idx = pool.size();
    pool.push_back(Quad());
    pool.push_back(Quad());
    pool.push_back(Quad());
    pool.push_back(Quad());
    return idx;*/
}

//Quad* Quad::get_quad(uint32_t idx) {
//    return &pool[idx];
//}