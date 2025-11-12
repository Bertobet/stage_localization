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

    std::vector<int> my_anchor_keys;
    if (node.storage(is_anchor{}))
        my_anchor_keys = { id };

    // Creazione hop map

    auto hop_map_all = spawn(CALL, [&](int nodeid){
        using fcpp::coordination::abf_hops;
        bool is_source = (node.uid == nodeid);
        hops_t d = abf_hops(CALL, is_source);
        return make_tuple(static_cast<int>(d), true);
    }, my_anchor_keys);   
    node.storage(hop_map{}) = hop_map_all;

    std::stringstream ss;
    if (hop_map_all.empty()) {
        ss << "(vuota)";
    } else {
        for (auto const& [key, value] : hop_map_all) 
            ss << "(id: " << key << " hop: " << value << " )";
    }

    // trasmissione posizione ancore
    auto information_anchor_map = spawn(CALL, [&](int nodeid){
        using fcpp::coordination::abf_hops;
        bool is_source = (node.uid == nodeid);
        auto distance = abf_hops(CALL, is_source);
        auto xs = broadcast(CALL, distance, x);
        auto ys = broadcast(CALL, distance, y);
        return make_tuple(make_tuple(xs, ys), true);
    }, my_anchor_keys);

    for (auto const& [id, t] : information_anchor_map) {
        node.storage(anchor_x_map{})[id] = get<0>(t);
        node.storage(anchor_y_map{})[id] = get<1>(t);
    }

    // Creazione map distance

    auto& xs_map = node.storage(anchor_x_map{});
    auto& ys_map = node.storage(anchor_y_map{});

    std::unordered_map<int,double, fcpp::common::hash<int>> distance_map_all;

    for (auto const& [id_anchor, x_anchor] : xs_map) {
        double y_anchor = ys_map.at(id_anchor);

        double dx = x_anchor - x;
        double dy = y_anchor - y;
        double dist = std::sqrt(dx*dx + dy*dy);

        distance_map_all[id_anchor] = dist;
    }

    node.storage(anchor_distance_map{}) = distance_map_all;

    std::stringstream supp2;
    if (node.storage(anchor_distance_map{}).empty()) {
        supp2 << "(vuota)";
    } else {
        for (auto const& [key, value] : node.storage(anchor_distance_map{})) 
            supp2 << "(id: " << key << " dista: " << value << " )";
    }


    // Calcolo correction

    int correction = 0;
    int hop= 0, distance = 0;
    if (node.storage(is_anchor{})){
        for (auto const& [key, value] : node.storage(anchor_distance_map{})){
            auto it = node.storage(hop_map{}).find(key);
            distance += value;
            hop += it->second;
        }   
        if (distance != 0 && hop != 0)
            node.storage(correction_anchor{}) = distance/hop;
    }

    // Trasmissione correction

    if (node.storage(is_anchor{}))
        node.storage(anchor_correction_map{})[id] = node.storage(correction_anchor{});
        
    auto correction_map_all = spawn(CALL, [&](int nodeid){
            using fcpp::coordination::abf_hops;
            bool is_source = (node.uid == nodeid);
            auto distance = abf_hops(CALL, is_source);
            int d = broadcast(CALL, distance, node.storage(correction_anchor{}));
            return make_tuple(d, true);
        }, my_anchor_keys);

    node.storage(anchor_correction_map{}) = correction_map_all;

    /*std::stringstream supp2;
    if (node.storage(anchor_y_map{}).empty()) {
        supp2 << "(vuota)";
    } else {
        for (auto const& [key, value] : node.storage(anchor_y_map{})) 
            supp2 << "(id: " << key << " y: " << value << " )";
    }*/

    //Calcolo distanza nodo-ancora
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

    //trilaterazione

    auto& dist_map = node.storage(distance_nodo_ancora_map{});
    auto& x_map = node.storage(anchor_x_map{});
    auto& y_map = node.storage(anchor_y_map{});

    std::vector<int> anchor_ids;
    for (auto const& [id, _] : dist_map)
        anchor_ids.push_back(id);

    int a1 = -1, a2 = -1, a3 = -1;

    if (anchor_ids.size() >= 3) {
    // Prendo la prima ancora come riferimento
    int a_ref = anchor_ids[0];
    double x1 = x_map[a_ref];
    double y1 = y_map[a_ref];
    double r1 = dist_map[a_ref];

    // Inizializzo i coefficienti della matrice normale
    double ATA11 = 0, ATA12 = 0, ATA22 = 0;
    double ATb1 = 0, ATb2 = 0;

    // Costruisci e accumula (AᵀA) e (Aᵀb)
    for (size_t i = 1; i < anchor_ids.size(); ++i) {
        int ai = anchor_ids[i];
        double xi = x_map[ai];
        double yi = y_map[ai];
        double ri = dist_map[ai];

        double Ai1 = 2 * (xi - x1);
        double Ai2 = 2 * (yi - y1);
        double bi  = r1*r1 - ri*ri + xi*xi - x1*x1 + yi*yi - y1*y1;

        ATA11 += Ai1 * Ai1;
        ATA12 += Ai1 * Ai2;
        ATA22 += Ai2 * Ai2;
        ATb1  += Ai1 * bi;
        ATb2  += Ai2 * bi;
    }

    // Determinante del sistema normale
    double det = ATA11 * ATA22 - ATA12 * ATA12;

    if (fabs(det) > 1e-12) {
        double x_est = (ATb1 * ATA22 - ATA12 * ATb2) / det;
        double y_est = (ATA11 * ATb2 - ATA12 * ATb1) / det;

        node.storage(x_stimato{}) = x_est;
        node.storage(y_stimato{}) = y_est;
    }
}


    // usage of node storage
    node.storage(node_size{})  = 10;
    node.storage(node_shape{}) = shape::sphere;
    //node.storage(spr{}) = supp.str();
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
    anchor_distance_map,    std::unordered_map<int,double, fcpp::common::hash<int>>,
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
