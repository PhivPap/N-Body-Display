#ifndef QUAD
#define QUAD

#include "util.h"
#include "body.h"
#include <vector>
#include <inttypes.h>

// The main quadtree class
class Quad {
private:
    void compute_bhtree_recursive(uint8_t depth);

public:
    /*static std::vector<Quad> pool;*/
    static Quad* pool;
    static uint32_t pool_size;
    static uint32_t pool_idx;

    Area area; 
    uint32_t body_count;
    double mass, diag_len_2;
    Point center_of_mass;
    std::vector<Body*> contained_bodies;
    
    // Children of this tree
    //uint32_t top_left_idx;
    Quad* top_left_quad;
    Quad* top_right_quad;
    Quad* bot_left_quad;
    Quad* bot_right_quad;

    Quad(void);
    Quad(Body* bodies, uint32_t body_count, const Area& area);
    ~Quad();
    void insert_body(Body* body);
    void set_area(const Area& area);
    static void set_pool(uint32_t init_size);
    static void reset_pool(void);
    static uint32_t pool_get_idx(void);
    /*static Quad* get_full_quad(uint32_t idx);*/
};


#endif