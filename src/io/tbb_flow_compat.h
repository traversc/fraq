#ifndef _FRAQ_TBB_FLOW_COMPAT_H
#define _FRAQ_TBB_FLOW_COMPAT_H

#include <utility>
#include <tbb/flow_graph.h>

// TBB 2019 uses source_node/decrement and surfaces TBB_INTERFACE_VERSION via
// tbb_stddef.h. oneTBB uses input_node/decrementer(), but <tbb/flow_graph.h>
// wrappers do not always surface the version macro. Prefer the interface
// version when the active flow-graph include exposes it, otherwise fall back
// to the old-header guard instead of probing unrelated oneapi/ headers that
// may come from a different package on the include path. Callers can override
// this probe if they know better.
#ifndef FRAQ_TBB_FLOW_USES_INPUT_NODE
#if defined(TBB_INTERFACE_VERSION)
#if TBB_INTERFACE_VERSION >= 12000
#define FRAQ_TBB_FLOW_USES_INPUT_NODE 1
#else
#define FRAQ_TBB_FLOW_USES_INPUT_NODE 0
#endif
#elif defined(__TBB_tbb_stddef_H)
#define FRAQ_TBB_FLOW_USES_INPUT_NODE 0
#else
#define FRAQ_TBB_FLOW_USES_INPUT_NODE 1
#endif
#endif

namespace fraq_internal {
namespace tbb_compat {

// Normalize the flow-graph API split between classic TBB (interface 11)
// and oneTBB (interface 12+).
//
// Classic TBB:
//   tbb::flow::source_node<T>
//   limiter_node.decrement
//   body signature: bool(T & out)
//
// oneTBB:
//   tbb::flow::input_node<T>
//   limiter_node.decrementer()
//   body signature: T(tbb::flow_control & fc)
//
// Normalized contract used below:
//   bool(T & out)
//
// Return true and fill `out` to emit one message.
// Return false to stop the node.
#if FRAQ_TBB_FLOW_USES_INPUT_NODE
template <class T>
struct source_node : tbb::flow::input_node<T> {
    using base_type = tbb::flow::input_node<T>;

    template <class Body>
    source_node(tbb::flow::graph &graph, Body &&body) :
        base_type(
            graph,
            [body = std::forward<Body>(body)](tbb::flow_control &fc) mutable -> T {
                T out{};
                if (!body(out)) {
                    fc.stop();
                }
                return out;
            }
        ) {}
};

template <class LimiterNode>
inline auto &decrementer(LimiterNode &node) {
    return node.decrementer();
}
#else
template <class T>
struct source_node : tbb::flow::source_node<T> {
    using base_type = tbb::flow::source_node<T>;

    template <class Body>
    source_node(tbb::flow::graph &graph, Body &&body) :
        base_type(graph, std::forward<Body>(body), false) {}
};

template <class LimiterNode>
inline auto &decrementer(LimiterNode &node) {
    return node.decrement;
}
#endif

} // namespace tbb_compat
} // namespace fraq_internal

#endif
