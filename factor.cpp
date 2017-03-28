/**
 * A simple driver program for factoring.  Doesn't really have error checking.
**/

#include <iostream>
#include <algorithm>
#include <vector>
#include <set>
#include <queue>
#include <math.h>

#include "bigint.h"

using namespace std;

void print_usage() {
  cout << "./factor [OPTIONS] [NUMBER]" << endl;
  cout << endl;
  cout << "Outputs a sorted list of the prime factors of NUMBER" << endl;
  cout << endl;
  cout << "[OPTIONS] may include:" << endl;
  cout << "  -r b - factor a random number with b bits.  [NUMBER] will be ignored if present" << endl;
  cout << "  -R b - factor a random semiprime with b bits.  [NUMBER] will be ignored" << endl;
  cout << "  -v   - give verbose output" << endl;
  cout << "  -t   - output timing information" << endl;
}

int main(int argc, const char * argv []) {
  srand(time(NULL));
  bigint n;

  int haven = 0;
  bool verbose = false;
  bool timing = false;
  for(int i = 1; i < argc; i++) {
    if(argv[i][0] == '-') {
      for(int j = 1; argv[i][j]; j++) {
        if(argv[i][j] == 'v') {
          verbose = true;
        } else if(argv[i][j] == 't') {
          timing = true;
        } else if(argv[i][j] == 'r') {
          if(haven) {
            cerr << "Two or more r/R flags cannot be present" << endl;
            return -1;
          }
          haven = 1;
        } else if(argv[i][j] == 'R') {
          if(haven) {
            cerr << "Two or more r/R flags cannot be present" << endl;
            return -1;
          }
          haven = 2;
        }
      }
      if(haven > 0) {
        if(i + 1 == argc) {
          cerr << "Expected number after r/R flag" << endl;
          return -1;
        }
        int bits = atoi(argv[++i]);
        if(haven == 1) {
          n = bigint::random(bits);
        } else if(haven == 2) {
          n = bigint::random_prime(bits / 2) *
              bigint::random_prime((bits + 1) / 2);
        }
        haven = -1;
      }
    }
  }
  if(!haven) {
    n = bigint(argv[argc - 1]);
  } else if(verbose) {
    cout << "Factoring " << n << endl;
  }

  clock_t start_time = clock();
  vector<bigint> f = n.factor(verbose);
  if(timing) {
    cout << "Factoring took " << 1.0 * (clock() - start_time) / CLOCKS_PER_SEC
         << " seconds" << endl;
  }

  for(int i = 0; i < f.size(); i++) {
    cout << f[i] << endl;
  }
  return 0;
}
