// Copyright © 2021 Giorgio Audrito. All Rights Reserved.

/**
 * @file exercises.cpp
 * @brief Quick-start aggregate computing exercises.
 */

// [INTRODUCTION]
//! Importing the FCPP library.
#include "lib/fcpp.hpp"

/**
 * @brief Namespace containing all the objects in the FCPP library.
 */
namespace fcpp {

//! @brief Dummy ordering between positions (allows positions to be used as secondary keys in ordered tuples).
template <size_t n>
bool operator<(vec<n> const&, vec<n> const&) {
    return false;
}

//! @brief Namespace containing the libraries of coordination routines.
namespace coordination {

//! @brief Tags used in the node storage.
namespace tags {
    //! @brief Color of the current node.
    struct node_color {};
    //! @brief Size of the current node.
    struct node_size {};
    //! @brief Shape of the current node.
    struct node_shape {};
    //! @brief is anchor or not
    struct is_anchor {};
    //! @brief hop map
    struct hop_map {};
    //! @brief support string
    struct spr {};
    //! @brief anchor distance map
    struct anchor_distance_map {};
    //! @brief correction for anchor
    struct correction_anchor {};
    //! @brief flag booleano
    struct flag_correction {};
    //! @brief map correction
    struct anchor_correction_map {};
    //! @brief map x
    struct anchor_x_map {};
    //! @brief map y
    struct anchor_y_map {};
    //! @brief distance nodo-ancora
    struct distance_nodo_ancora_map {};
    //! @brief x stimato
    struct x_stimato {};
    //! @brief y stimato
    struct y_stimato {};
    //! @brief support string
    struct sprr {};
}

//! @brief The maximum communication range between nodes.
constexpr size_t communication_range = 100;

// [AGGREGATE PROGRAM]

/*
 * @brief example function for checking a property. 
 * Sample property: you (the current device) have not been at disrisk for a couple of rounds.
 */
FUN bool recent_dis_monitor(ARGS, bool disrisk) { CODE
    using namespace logic;
    bool prev_disrisk = Y(CALL, disrisk);
    return !disrisk & !prev_disrisk;
}
FUN_EXPORT monitor_t = export_list<past_ctl_t, slcs_t>;

// @brief Main function.
MAIN() {
    // import tag names in the local scope.
    using namespace tags;

    int id = node.uid;
    double side = 500.0;  
    double step = 100.0;

    int anchors_per_side = static_cast<int>(side / step);
    int total_anchors = anchors_per_side * 4; 

    double x = 0, y = 0;

    if (id < total_anchors) {
        int pos = id;

        if (pos <= anchors_per_side) {
            x = pos * step;
            y = side;
        }
        else if (pos <= anchors_per_side * 2) {
            x = side;
            y = side - (pos - anchors_per_side) * step;
        }
        else if (pos <= anchors_per_side * 3) {
            x = side - (pos - anchors_per_side * 2) * step;
            y = 0;
        }
        else {
            x = 0;
            y = (pos - anchors_per_side * 3) * step;
        }

        node.position() = make_vec(x, y);
        node.storage(is_anchor{}) = true;
        node.storage(node_color{}) = color(RED);
    } else {
        node.storage(is_anchor{}) = false;
        node.storage(node_color{}) = color(GREEN);
    }

    // Creazione hop map

    std::vector<int> my_keys = { id };
    auto hop_map_all = spawn(CALL, [&](int nodeid){
        using fcpp::coordination::abf_hops;
        bool is_source = (node.uid == nodeid);
        hops_t d = abf_hops(CALL, is_source);
        return make_tuple(static_cast<int>(d), true);
    }, my_keys);   
    node.storage(hop_map{}) = hop_map_all;

    std::stringstream ss;
    if (hop_map_all.empty()) {
        ss << "(vuota)";
    } else {
        for (auto const& [key, value] : hop_map_all) 
            ss << "(id: " << key << " hop: " << value << " )";
    }

    // Creazione map distance

    std::vector<int> my_anchor_keys;
    if (node.storage(is_anchor{}))
        my_anchor_keys = { id };

    auto distance_map_all = spawn(CALL, [&](int nodeid){
        using fcpp::coordination::abf_hops;
        bool is_source = (node.uid == nodeid);
        int d = abf_distance(CALL, is_source);
        return make_tuple(d, true);
    }, my_anchor_keys);

    node.storage(anchor_distance_map{}) = distance_map_all;

    std::stringstream supp2;
    if (distance_map_all.empty()) {
        supp2 << "(vuota)";
    } else {
        for (auto const& [key, value] : distance_map_all) 
            supp2 << "(id: " << key << " dista: " << value << " )";
    }

    // Calcolo correction


    int correction = 0;
    int hop= 0, distance = 0;
    if (node.storage(is_anchor{}) && !node.storage(flag_correction{}) && node.current_time() > 10){
        for (auto const& [key, value] : node.storage(anchor_distance_map{})){
            auto it = node.storage(hop_map{}).find(key);
            distance += value;
            hop += it->second;
        }   
        node.storage(flag_correction{}) = true;
        node.storage(correction_anchor{}) = distance/hop;
    }

    // Trasmissione correction

    if ((node.storage(is_anchor{}) && node.storage(flag_correction{})) || node.current_time() > 11){
        node.storage(anchor_correction_map{})[id] = node.storage(correction_anchor{});
        

        auto correction_map_all = spawn(CALL, [&](int nodeid){
            using fcpp::coordination::abf_hops;
            bool is_source = (node.uid == nodeid);
            auto distance = abf_hops(CALL, is_source);
            int d = broadcast(CALL, distance, node.storage(correction_anchor{}));
            return make_tuple(d, true);
        }, my_anchor_keys);

        node.storage(anchor_correction_map{}) = correction_map_all;
    }
    


    /*std::stringstream supp;
    if (node.storage(anchor_correction_map{}).empty()) {
        supp << "(vuota)";
    } else {
        for (auto const& [key, value] : node.storage(anchor_correction_map{})) 
            supp << "(id: " << key << " correction: " << value << " )";
    }*/

    // Trasmissione position ancore

    if (node.storage(is_anchor{}))
        node.storage(anchor_x_map{})[id] = x;

    auto x_map_all = spawn(CALL, [&](int nodeid){
        using fcpp::coordination::abf_hops;
        bool is_source = (node.uid == nodeid);
        auto distance = abf_hops(CALL, is_source);
        int d = broadcast(CALL, distance, x);
        return make_tuple(d, true);
    }, my_anchor_keys);

    node.storage(anchor_x_map{}) = x_map_all;

    if (node.storage(is_anchor{}))
        node.storage(anchor_y_map{})[id] = y;

    auto y_map_all = spawn(CALL, [&](int nodeid){
        using fcpp::coordination::abf_hops;
        bool is_source = (node.uid == nodeid);
        auto distance = abf_hops(CALL, is_source);
        int d = broadcast(CALL, distance, y);
        return make_tuple(d, true);
    }, my_anchor_keys);

    node.storage(anchor_y_map{}) = y_map_all;

    //Calcolo distanza nodo-ancora
    if (node.current_time() > 13){
        for (int i = 0; i < total_anchors; i++){
            auto& hop_map_ref = node.storage(hop_map{});
            auto& corr_map_ref = node.storage(anchor_correction_map{});
            auto& dist_map_ref = node.storage(distance_nodo_ancora_map{});

            auto hop_it = hop_map_ref.find(i);
            auto corr_it = corr_map_ref.find(i);

            //Controlla che entrambe le chiavi esistano
            if (hop_it != hop_map_ref.end() && corr_it != corr_map_ref.end()) {
                dist_map_ref[i] = hop_it->second * corr_it->second;
            }
        }
    }
      
    std::stringstream supp;
    if (node.storage(distance_nodo_ancora_map{}).empty()) {
        supp << "(vuota)";
    } else {
        for (auto const& [key, value] : node.storage(distance_nodo_ancora_map{})) 
            supp << "(id: " << key << " correction: " << value << " )";
    }

    //trilaterazione

    auto& dist_map = node.storage(distance_nodo_ancora_map{});
    auto& x_map = node.storage(anchor_x_map{});
    auto& y_map = node.storage(anchor_y_map{});

    std::vector<int> anchor_ids;
    for (auto const& [id, _] : dist_map)
        anchor_ids.push_back(id);

    int a1 = -1, a2 = -1, a3 = -1;
    bool found = false;

    // Cerca 3 ancore non allineate
    if (anchor_ids.size() >= 3) {
        for (size_t i = 0; i < anchor_ids.size() && !found; ++i) {
            for (size_t j = i + 1; j < anchor_ids.size() && !found; ++j) {
                for (size_t k = j + 1; k < anchor_ids.size() && !found; ++k) {
                    double x1 = x_map[anchor_ids[i]];
                    double y1 = y_map[anchor_ids[i]];
                    double x2 = x_map[anchor_ids[j]];
                    double y2 = y_map[anchor_ids[j]];
                    double x3 = x_map[anchor_ids[k]];
                    double y3 = y_map[anchor_ids[k]];

                    // Calcolo area per verificare che non siano allineate
                    double area = fabs((x2 - x1)*(y3 - y1) - (y2 - y1)*(x3 - x1)) / 2.0;
                    if (area > 1e-6) {
                        a1 = anchor_ids[i];
                        a2 = anchor_ids[j];
                        a3 = anchor_ids[k];
                        found = true;
                    }
                }
            }
        }
    }

    if (found) {
        double x1 = x_map[a1], y1 = y_map[a1], r1 = dist_map[a1];
        double x2 = x_map[a2], y2 = y_map[a2], r2 = dist_map[a2];
        double x3 = x_map[a3], y3 = y_map[a3], r3 = dist_map[a3];

        // --- Sistema lineare per la trilaterazione ---
        double A11 = 2*(x2 - x1);
        double A12 = 2*(y2 - y1);
        double A21 = 2*(x3 - x1);
        double A22 = 2*(y3 - y1);

        double b1 = r1*r1 - r2*r2 + x2*x2 - x1*x1 + y2*y2 - y1*y1;
        double b2 = r1*r1 - r3*r3 + x3*x3 - x1*x1 + y3*y3 - y1*y1;

        double det = A11*A22 - A12*A21;

        if (fabs(det) > 1e-12) {
            double x_est = ( b1*A22 - A12*b2 ) / det;
            double y_est = (-b1*A21 + A11*b2 ) / det;

            node.storage(x_stimato{}) = x_est;
            node.storage(y_stimato{}) = y_est;
        }
    }




    /*if (node.storage(distance_nodo_ancora_map{}).size() >= 3 && node.current_time() > 15){
        auto& dist_map = node.storage(distance_nodo_ancora_map{});
        auto& x_map = node.storage(anchor_x_map{});
        auto& y_map = node.storage(anchor_y_map{});

        auto it = dist_map.begin();
        int a1 = it->first; double x1 = x_map[a1]; double y1 = y_map[a1]; double r1 = it->second;
        ++it;
        int a2 = it->first; double x2 = x_map[a2]; double y2 = y_map[a2]; double r2 = it->second;
        ++it;
        int a3 = it->first; double x3 = x_map[a3]; double y3 = y_map[a3]; double r3 = it->second;

        double A11 = 2 * (x2 - x1);
        double A12 = 2 * (y2 - y1);
        double A21 = 2 * (x3 - x1);
        double A22 = 2 * (y3 - y1);

        double b1 = r1*r1 - r2*r2 + x2*x2 - x1*x1 + y2*y2 - y1*y1;
        double b2 = r1*r1 - r3*r3 + x3*x3 - x1*x1 + y3*y3 - y1*y1;

        double det = A11*A22 - A12*A21;

        if (fabs(det) > 1e-12) { // se le ancore non sono allineate
            double x_est = ( b1*A22 - A12*b2 ) / det;
            double y_est = (-b1*A21 + A11*b2 ) / det;
            node.storage(x_stimato{}) = x_est;
            node.storage(y_stimato{}) = y_est;
        } else {
            // Caso degenerato
            node.storage(x_stimato{}) = 0;
            node.storage(y_stimato{}) = 0;
        }

    }*/


    // usage of node storage
    node.storage(node_size{})  = 10;
    node.storage(node_shape{}) = shape::sphere;
    node.storage(spr{}) = supp.str();
     node.storage(sprr{}) = supp2.str();
    

}
//! @brief Export types used by the main function (update it when expanding the program).
FUN_EXPORT main_t = export_list<double, int, monitor_t/*, hop_from_node_t*/>;

} // namespace coordination

// [SYSTEM SETUP]

//! @brief Namespace for component options.
namespace option {

//! @brief Import tags to be used for component options.
using namespace component::tags;
//! @brief Import tags used by aggregate functions.
using namespace coordination::tags;

//! @brief Number of people in the area.
constexpr int node_num = 100;
//! @brief Dimensionality of the space.
constexpr size_t dim = 2;

//! @brief Description of the round schedule.
using round_s = sequence::periodic<
    distribution::interval_n<times_t, 0, 1>,    // uniform time in the [0,1] interval for start
    distribution::weibull_n<times_t, 10, 1, 10> // weibull-distributed time for interval (10/10=1 mean, 1/10=0.1 deviation)
>;
//! @brief The sequence of network snapshots (one every simulated second).
using log_s = sequence::periodic_n<1, 0, 1>;
//! @brief The sequence of node generation events (node_num devices all generated at time 0).
using spawn_s = sequence::multiple_n<node_num, 0>;
//! @brief The distribution of initial node positions (random in a 500x500 square).
using rectangle_d = distribution::rect_n<1, 0, 0, 500, 500>;
//! @brief The contents of the node storage as tags and associated types.
using store_t = tuple_store<
    node_color,             color,
    node_size,              double,
    node_shape,             shape,
    is_anchor,              bool,
    hop_map,                std::unordered_map<int,int, fcpp::common::hash<int>>,
    spr,                    std::string,
    anchor_distance_map,    std::unordered_map<int,int, fcpp::common::hash<int>>,
    correction_anchor,      int,
    flag_correction,        bool,
    anchor_correction_map,  std::unordered_map<int,int, fcpp::common::hash<int>>,
    anchor_x_map,           std::unordered_map<int,int, fcpp::common::hash<int>>,
    anchor_y_map,           std::unordered_map<int,int, fcpp::common::hash<int>>,
    distance_nodo_ancora_map, std::unordered_map<int,int, fcpp::common::hash<int>>,
    x_stimato,              double,
    y_stimato,              double,
    sprr,                   std::string
>;
//! @brief The tags and corresponding aggregators to be logged (change as needed).
using aggregator_t = aggregators<
    node_size, aggregator::mean<double>
>;

//! @brief The general simulation options.
DECLARE_OPTIONS(list,
    parallel<true>,      // multithreading enabled on node rounds
    synchronised<false>, // optimise for asynchronous networks
    program<coordination::main>,   // program to be run (refers to MAIN above)
    exports<coordination::main_t, std::unordered_set<int, fcpp::common::hash<int>>, std::unordered_map<int, int, fcpp::common::hash<int>>, tuple<int, int>>, // export type list (types used in messages)
    retain<metric::retain<2,1>>,   // messages are kept for 2 seconds before expiring
    round_schedule<round_s>, // the sequence generator for round events on nodes
    log_schedule<log_s>,     // the sequence generator for log events on the network
    spawn_schedule<spawn_s>, // the sequence generator of node creation events on the network
    store_t,       // the contents of the node storage
    aggregator_t,  // the tags and corresponding aggregators to be logged
    init<
        x,      rectangle_d // initialise position randomly in a rectangle for new nodes
    >,
    dimension<dim>, // dimensionality of the space
    connector<connect::fixed<100, 1, dim>>, // connection allowed within a fixed comm range
    shape_tag<node_shape>, // the shape of a node is read from this tag in the store
    size_tag<node_size>,   // the size  of a node is read from this tag in the store
    color_tag<node_color>  // the color of a node is read from this tag in the store
);

} // namespace option

} // namespace fcpp


//! @brief The main function.
int main() {
    using namespace fcpp;

    //! @brief The network object type (interactive simulator with given options).
    using net_t = component::interactive_simulator<option::list>::net;
    //! @brief The initialisation values (simulation name).
    auto init_v = common::make_tagged_tuple<option::name>("Exercises");
    //! @brief Construct the network object.
    net_t network{init_v};
    //! @brief Run the simulation until exit.
    network.run();
    return 0;
}
