/**
 * @file boyermoore.cpp
 * @author Jiawei Huang
 * @brief implementation of Boyer_Moore algorithm
 * @version 0.1
 * @date 2020-02-07
 *
 * @copyright Copyright (c) 2020
 *
 */

#include <algorithm>
#include <cassert>
#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <vector>

size_t dna_encoding[] = {
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0,
    4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
};

class Boyer_Moore {
public:
  /**
   * @brief:    map print function
   *
   * @param mymap:      map to be printed, values are of type size_t
   */
  void print_map(std::map<char, size_t> &mymap) {
    for (std::map<char, size_t>::iterator it = mymap.begin(); it != mymap.end();
         ++it)
      std::cout << it->first << " => " << it->second << '\n';
  }

  /**
   * @brief:     override for print map function
   *
   * @param mymap:      map to be printed, values are vectors
   */
  void print_map(std::map<char, std::vector<size_t>> &mymap) {
    for (std::map<char, std::vector<size_t>>::iterator it = mymap.begin();
         it != mymap.end(); ++it) {
      std::cout << it->first << " : ";
      print_arr(it->second);
    }
  }

  /**
   * @brief     print vectors
   *
   * @param arr     vector to be printed
   */
  void print_arr(std::vector<size_t> &arr) {
    for (auto &i : arr)
      std::cout << i << " ";
    std::cout << '\n';
  }

  void print_mismatch(const std::string &pat, const std::string &text,
                      const size_t &i, const size_t &j, const size_t &shift) {
    size_t n = text.size(), m = pat.size();
    std::cout << text.substr(0, i + j) << " | " << text.substr(i + j, 1)
              << " | " << text.substr(i + j + 1, n - i - j) << '\n';

    std::cout << std::string(i, ' ') << pat.substr(0, j) << " | "
              << pat.substr(j, 1) << " | " << pat.substr(j + 1, m - j) << '\n';

    std::cout << "after shifting " << shift << " space text starts matching at "
              << text.substr(0, i + shift) << " | "
              << text.substr(i + shift, n - i - shift) << '\n';
  }
  //***********************************************
  //     1. bad character heuristic
  //***********************************************
  /**
   * @brief     return skips of bad character rule when mismatch happend
   *
   * @param tab   pattern shift table
   * @param text  mismatch character happend in text
   * @param j     mismatch index in pattern
   * @return size_t   shift space
   */
  size_t bad_character_mismatch(std::vector<std::vector<size_t>> &tab,
                                const char &text, const size_t &j) {
    size_t shift = 1;
    if (dna_encoding[text] > 3)
      shift = j + 1;
    else {
      size_t i = dna_encoding[text];
      shift = j - tab[j][i] + 1;
    }
    return shift;
  }

  /**
   * @brief     Given pattern string and list with ordered alphabet characters,
   create and return a dense bad character table.  Table is indexed by offset
                then by character.
   *
   * @param p     pattern
   * @return std::vector<std::vector<size_t>> bad character table
   */
  std::vector<std::vector<size_t>> get_bad_char_table(const std::string &p) {
    std::vector<size_t> index(p.size(), 0);
    std::vector<std::vector<size_t>> table;

    for (size_t i = 0; i < p.size(); ++i) {
      char c = p[i];
      table.push_back(index);
      index[dna_encoding[c]] = i + 1;
    }

    return table;
  }

  //*******************************************
  //     2. goof suffix rule
  //*******************************************
  /**
   * @brief     naive match function, one string with two different starting
   * positions
   *
   * @param s       matching string
   * @param q       first position index
   * @param i       second position index
   * @return size_t         matching end position of former one(q)
   */
  static size_t match(const std::string &s, size_t q, size_t i) {
    while (i < s.length() && s[q] == s[i]) {
      ++q;
      ++i;
    }
    return q;
  }

  /**
   * @brief       Z value computation algorithm (Andrew's version)
   *
   * @param s     string
   * @return std::vector<size_t>      array storing z values of string s
   */
  std::vector<size_t> z_array(const std::string &s) {
    std::vector<size_t> Z(s.length());

    size_t l = 0, r = 0;
    for (size_t k = 1; k < s.length(); ++k) {
      if (k >= r) { // Case 1: full comparison
        Z[k] = match(s, 0, k);
        r = k + Z[k];
        l = k;
      } else { // Case 2: (we are inside a Z-box)
        const size_t k_prime = k - l;
        const size_t beta_len = r - k;
        if (Z[k_prime] < beta_len) { // Case 2a: stay inside Z-box
          Z[k] = Z[k_prime];
        } else { // Case 2b: need to match outside the Z-box
          Z[k] = match(s, beta_len, r);
          r = k + Z[k];
          l = k;
        }
      }
    }
    return Z;
  }

  /**
   * @brief         To calculate the N array for good suffix rule,
   *                N[j] = Z[n-j+1], Z is the z-value list of reverse s
   *
   * @param s       pattern string
   * @return std::vector<size_t>        N array
   */
  std::vector<size_t> n_array(const std::string &s) {
    std::string sr(s);
    // reverse s to get z value of reversed pattern
    std::reverse(sr.begin(), sr.end());

    std::vector<size_t> z = z_array(sr);
    // reverse again to get N value. Z[i] = N[n-j+1], e.g. N[1] = Z[n], N[2] =
    // Z[n-1]
    std::reverse(z.begin(), z.end());

    return z;
  }

  /**
   * @brief     This function is to calculate the L' array for good suffix rule
   *            L'[i] = largest index j less than n such that N[j] = |P[i:n]|
   *
   * @param p       pattern string
   * @param n       N array got from reverse function n_array
   * @return std::vector<size_t>        L' array
   */
  std::vector<size_t> big_l_prime_array(const std::string &p,
                                        std::vector<size_t> &n) {
    std::vector<size_t> lp(p.size(), 0);
    for (size_t j = 0; j < p.size() - 1; ++j) {
      // cause index begin at 0, so the (p.size() - n[j] + 1)th elem's index
      // should be p.size() - n[j]
      size_t i = p.size() - n[j];
      if (i < p.size())
        // j is the index of n(counting starts at 0), but we need to store the
        // value of lp, so we use j + 1
        lp[i] = j + 1;
    }
    return lp;
  }

  /**
   * @brief     calculate the L array for good suffix rule
   *            L[i] = largest index j less than n such that N[j] >= |P[i:n]|
   *
   * @param p       pattern
   * @param lp      L' array get from big_l_prime_array
   * @return std::vector<size_t>        L array
   */
  std::vector<size_t> big_l_array(std::string p, std::vector<size_t> lp) {
    std::vector<size_t> l(p.size(), 0);
    l[1] = lp[1];
    for (size_t i = 2; i < p.size(); ++i)
      l[i] = std::max(l[i - 1], lp[i]);

    return l;
  }

  /**
   * @brief         calculate the l' array for good suffix rule
   *                l'[i] = largest j <= |P[i:n]| = n-i+1, such that N[j] = j
   * @param n       N array got from reverse z algorithm
   * @return std::vector<size_t>        l' array
   */
  std::vector<size_t> small_l_prime_array(std::vector<size_t> &n) {
    std::vector<size_t> small_lp(n.size(), 0);

    // first get the rightmost N[j] = j, then i = n-j+1(cause index starts at 0,
    // so here is N[j] = j+1, i = n-j-1)
    for (size_t i = 0; i < n.size(); ++i) {
      if (n[i] == i + 1) // prefix matching a suffix
        small_lp[n.size() - i - 1] = i + 1;
    }
    // then for any i < j, l'[i] = l'[j]
    for (int i = n.size() - 2; i >= 0; --i) {
      if (small_lp[i] == 0) // smear them out to the left
        small_lp[i] = small_lp[i + 1];
    }

    return small_lp;
  }

  /**
   * @brief         Given a mismatch at offset i, and given L/L' and l' arrays,
   *
   * @param i       mismatch index
   * @param big_l_prime      L' array get from big_l_prime_array
   * @param small_l_prime       l' array get from small_l_prime_array
   * @return size_t        amount to shift as determined by good suffix rule.
   */
  size_t good_suffix_mismatch(int i, const std::vector<size_t> &big_l_prime,
                              const std::vector<size_t> &small_l_prime) {
    size_t length = big_l_prime.size();
    assert(i < length);

    if (i == length - 1)
      return 1;

    i++; // i points to the leftmost matching position of P
    if (big_l_prime[i] > 0)
      return (length - big_l_prime[i]);

    return (length - small_l_prime[i]);
  }

  /**
   * @brief       boyer moore based on Gusfield
   *
   * @param p         pattern
   * @param t         text
   * @return size_t       number of matching occurrence
   */
  void boyer_moore(const std::string &p, const std::string &t) {
    size_t i = 0, match = 0, comparision = 0;
    std::vector<std::vector<size_t>> tab = get_bad_char_table(p);

    std::vector<size_t> narray = n_array(p);
    std::vector<size_t> biglprime = big_l_prime_array(p, narray);
    std::vector<size_t> bigl = big_l_array(p, biglprime);
    std::vector<size_t> smalllprime = small_l_prime_array(narray);

    size_t m = p.length(), n = t.length();
    auto start = std::chrono::steady_clock::now();
    while (i < n - m + 1) {
      size_t shift = 1;
      bool mismatch = false;
      size_t j = m;
      while (j--) {
        if (p[j] != t[i + j]) {
          size_t skip_bc = bad_character_mismatch(tab, t[i + j], j);
          size_t skip_gs = good_suffix_mismatch(j, biglprime, smalllprime);
          shift = std::max(skip_bc, skip_gs);
          // print mismatch position
          // std::cout << "\nmismatch happend, bad character shift: " << skip_bc
          //           << " , good suffix shift: " << skip_gs << '\n';
          // print_mismatch(p, t, i, j, shift);
          mismatch = true;
          comparision += m - j;
          break;
        }
      }
      if (mismatch == false) {
        match++;
        // print matching position
        // std::cout << "\npattern matched at index " << i << " with text "
        //           << t.substr(0, i) << " " << t.substr(i, m) << " "
        //           << t.substr(i + m, n - i - m) << '\n';
        size_t skip_gs = smalllprime.size() - smalllprime[1];
        shift = std::max(shift, skip_gs);
        comparision += m;
      }
      i += shift;
    }
    auto end = std::chrono::steady_clock::now();
    auto diff = end - start;

    std::cout << "total number of matching: " << match << '\n';
  }
};

/**
 * @brief       reads a FASTA format file line-by-line, skipping the "name"
 * lines
 *
 * @param fasta_filename
 * @param T     string of text used in boyer moore
 */
static void read_fasta(const std::string &fasta_filename, std::string &T) {
  std::ifstream in(fasta_filename);
  if (!in)
    throw std::runtime_error("problem with file: " + fasta_filename);

  const std::streampos begin_pos = in.tellg();
  in.seekg(0, std::ios_base::end); // move to the end
  const size_t filesize = in.tellg() - begin_pos;
  in.seekg(0, std::ios_base::beg); // and now move to back to the beginning

  T.clear();           // make sure we start fresh
  T.reserve(filesize); // reserve enough total space
  std::string line;
  while (getline(in, line))
    if (line[0] != '>')
      T += line;
}

int main(int argc, const char *const argv[]) {
  if (argc != 3) {
    std::cerr << "usage: " << argv[0] << " <PATTERN> <TEXT>" << std::endl;
    return EXIT_FAILURE;
  }

  const std::string P(argv[1]);
  std::string T;
  read_fasta(argv[2], T);
  // std::cout << "length: " << T.length();

  // make sure pattern not bigger than text
  assert(P.length() <= T.length());

  // initialize
  Boyer_Moore bm = Boyer_Moore();
  
  bm.boyer_moore(P, T);
}
