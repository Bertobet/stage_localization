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

// Funzione del processo
/*FUN fcpp::tuple<int,bool> hop_from_node(ARGS, message const& m){ CODE
    using fcpp::coordination::abf_hops;

    bool is_source = (node.uid == 1);

    hops_t d = abf_hops(CALL, is_source);

    return make_tuple(static_cast<int>(d), true);
}
FUN_EXPORT hop_from_node_t = export_list<hops_t>;*/

// @brief Main function.
MAIN() {
    // import tag names in the local scope.
    using namespace tags;

    int id = node.uid;

    bool anchor = (id < 8);
    node.storage(is_anchor{}) = anchor;
    node.storage(node_color{}) = anchor ? color(RED) : color(GREEN);

    // Nel caso volessi che le ancore abbiano una posizione specifica
    //if (id == 0) node.position() = make_vec(0,   0);

    std::vector<int> my_keys = { id };
    //fcpp::tuple<int,bool> (*pointer)(ARGS, int) = hop_from_node;
    //auto hop_map_all = spawn(CALL, hop_from_node, my_keys);

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
            ss << "(id: " << key << " dista: " << value << " ) ";
    }


    // usage of node physics
    //node.velocity() = -node.position()/communication_range;

    // usage of node storage
    node.storage(node_size{})  = 10;
    node.storage(node_shape{}) = shape::sphere;
    node.storage(spr{}) = ss.str();

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
    node_color, color,
    node_size,  double,
    node_shape, shape,
    is_anchor,  bool,
    hop_map,    std::unordered_map<int,int, fcpp::common::hash<int>>,
    spr,        std::string
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
    exports<coordination::main_t, std::unordered_set<int, fcpp::common::hash<int>>>, // export type list (types used in messages)
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
