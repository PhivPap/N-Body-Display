#ifndef QUAD
#define QUAD

#include "MyTypes.h"
#include "util.h"
#include "body.h"
#include <vector>
#include <inttypes.h>

// The main quadtree class
class Quad {
private:
    void compute_bhtree_recursive(u8 depth);

public:
    /*static std::vector<Quad> pool;*/
    static Quad* pool;
    static u32 pool_size;
    static u32 pool_idx;

    Area area; 
    u32 body_count;
    f64 mass, diag_len_2;
    Point center_of_mass;
    std::vector<Body*> contained_bodies;
    
    // Children of this tree
    //u32 top_left_idx;
    Quad* top_left_quad;
    Quad* top_right_quad;
    Quad* bot_left_quad;
    Quad* bot_right_quad;

    Quad(void);
    Quad(Body* bodies, u32 body_count, const Area& area);
    ~Quad();
    void insert_body(Body* body);
    void set_area(const Area& area);
    static void set_pool(u32 init_size);
    static void reset_pool(void);
    static u32 pool_get_idx(void);
    /*static Quad* get_full_quad(u32 idx);*/
};


#endif