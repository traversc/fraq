#ifndef _FRAQ_NTHREADS_GUARD_H_
#define _FRAQ_NTHREADS_GUARD_H_

#if !defined(_WIN32)
#include <sys/types.h>
#include <unistd.h>
#endif

namespace fraq_internal {

#if defined(_WIN32)

inline void loaded_in_fork_child_internal(bool loaded_in_fork_child) {
  (void)loaded_in_fork_child;
}

inline bool is_forked_child_internal() {
  return false;
}

#else

// inline, not static: static would give every translation unit including this
// header its own copy, so a flag set through one TU would be invisible to the
// inline readers compiled into another.
inline const pid_t fraq_load_pid_internal = getpid();
inline bool fraq_loaded_in_fork_child_state_internal = false;

inline void loaded_in_fork_child_internal(bool loaded_in_fork_child) {
  fraq_loaded_in_fork_child_state_internal = loaded_in_fork_child;
}

inline bool is_forked_child_internal() {
  return fraq_loaded_in_fork_child_state_internal || getpid() != fraq_load_pid_internal;
}

#endif

inline int normalize_nthreads(int nthreads) {
  if (nthreads <= 1) {
    return 1;
  }
  if (is_forked_child_internal()) {
    return 1;
  }
  return nthreads;
}

} // namespace fraq_internal

#endif
