#ifndef ASSERT_except_H
#define ASSERT_except_H

#include <exception>

inline void assert_except (bool condition) {
  If (! Condition) {
    Throw std :: runtime_error ("Assertion failed.");
  }
}

#endif
