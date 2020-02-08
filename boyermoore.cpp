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
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <vector>

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
                      size_t &i, const size_t &j, const size_t &shift) {
    size_t n = text.size(), m = pat.size();
    std::cout << text.substr(0, i + j) << " | " << text.substr(i + j, 1)
              << " | " << text.substr(i + j + 1, n - i - j) << '\n';

    std::cout << std::string(i, ' ') << pat.substr(0, j) << " | "
              << pat.substr(j, 1) << " | " << pat.substr(j + 1, m - j) << '\n';

    i += shift;
    std::cout << "after shifting " << shift << " space text starts matching at "
              << text.substr(0, i) << " | " << text.substr(i, n - i) << '\n';
  }

  //     1. bad character heuristic

  /*
   * This function use the basic bad character rule
   * Just record the index of each character last showed.
   */

  /**
   * @brief   preprocess of bad character rule. Record the index of each
   *          character in pattern, if one character occurs more than twice,
   *          record the rightmost one
   *
   * @param pat   pattern
   * @param alphabet map to store index in pattern, for unique characters showed
   * in text while not in pattern, default 0
   */
  void processing(const std::string &pat, std::map<char, size_t> &alphabet) {
    for (size_t i = 0; i < pat.size(); ++i) {
      if (alphabet.find(pat[i]) != alphabet.end())
        alphabet[pat[i]] = i;
    }
    std::cout << "Bad character array list :\n";
    print_map(alphabet);
  }

  /**
   * @brief     preprocessing of extended bad character rule. Using a list to
   *            record all the occurrence of each character of the pattern
   *
   * @param pat     pattern
   * @param ex_list     list to store all the occurrence of characters
   */
  void extended_processing(const std::string &pat,
                           std::map<char, std::vector<size_t>> &ex_list) {
    // scan P from right to left(reversed order)
    size_t i = pat.size();
    while (i--) {
      if (ex_list.find(pat[i]) == ex_list.end())
        ex_list[pat[i]] = std::vector<size_t>(1, i);
      else
        ex_list[pat[i]].push_back(i);
    }

    std::cout << "character : index\n";
    print_map(ex_list);
  }

  /*
   * This function is to calculate the shift space for bad character rule
   * :param ex_list: extended list containing all occurence of each character in
   *                 pattern.
   * :param text:    text.
   * :param i:       start index of the text
   * : param j:      current index of the pattern
   */
  int get_bmshift(std::map<char, std::vector<size_t>> ex_list, std::string text,
                  size_t i, size_t j) {
    std::map<char, std::vector<size_t>>::iterator it =
        ex_list.find(text[i + j]);
    int shift = 0;
    if (it == ex_list.end())
      shift = j + 1;
    else {
      for (auto &k : it->second) {
        if (k < j) {
          shift = j - k;
          break;
        } else
          shift = j + 1;
      }
    }
    return shift;
  }

  void badchar(std::string text, std::string pat) {
    size_t n = text.size();
    size_t m = pat.size();
    std::map<char, size_t> alphabet;
    for (size_t i = 0; i < text.size(); ++i) {
      alphabet[text[i]] = 0;
    }
    processing(pat, alphabet);
    std::map<char, std::vector<size_t>> ex_list;
    extended_processing(pat, ex_list);

    int i = 0, j = m - 1;
    while (i <= n - m) {
      while (pat[j] != text[i + j]) {
        std::cout << "mismatch happend \n";
        int shift = get_bmshift(ex_list, text, i, j);
        // int shift = std::max(j - alphabet.at(text[i+j]), 1);
        i += shift;
        print_mismatch(pat, text, i, j, shift);
        j = m - 1;
      }
      j--;
      if (j < 0) {
        std::cout << "pattern matched at index " << i << " with text "
                  << text.substr(0, i) << " " << text.substr(i, m) << " "
                  << text.substr(i + m, n - i - m) << '\n';
        j = m - 1;
        i++;
      }
    }
  }

  //     2. goof suffix rule
  static size_t match(const std::string &s, size_t q, size_t i) {
    while (i < s.length() && s[q] == s[i]) {
      ++q;
      ++i;
    }
    return q;
  }

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
      // cout << k + 1 << "\t" << l + 1 << "\t" << r << "\t"
      //     << Z[k] << "\t" << the_case << endl;
    }
    return Z;
  }

  /*
   * This function is to calculate the N array for good suffix rule
   * N[j] = Z[n-j+1], Z is the z-value list of reverse s
   * :param s: pattern string
   */
  std::vector<size_t> n_array(const std::string &s) {
    std::string sr(s);
    std::reverse(sr.begin(), sr.end());

    std::vector<size_t> z = z_array(sr);
    std::reverse(z.begin(), z.end());

    return z;
  }

  /*
   * This function is to calculate the L' array for good suffix rule
   * L'[i] = largest index j less than n such that N[j] = |P[i:n]|
   * :param p: pattern string
   * :param n: N array got from reverse z algorithm
   */
  std::vector<size_t> big_l_prime_array(const std::string &p,
                                        std::vector<size_t> &n) {
    std::vector<size_t> lp(p.size(), 0);
    for (size_t j = 0; j < p.size() - 1; ++j) {
      // cause index begin at 0, so the (p.size() - n[j] + 1)th elem's index
      // should be p.size() - n[j]
      int i = p.size() - n[j];
      if (i < p.size())
        // j is the index of n(counting starts at 0), but we need to store the
        // value of lp, so we use j + 1
        lp[i] = j + 1;
    }
    return lp;
  }
  /*
   * This function is to calculate the L array for good suffix rule
   * L[i] = largest index j less than n such that N[j] >= |P[i:n]|
   * :param p: pattern string
   * :param lp: L' array
   */
  std::vector<size_t> big_l_array(std::string p, std::vector<size_t> lp) {
    std::vector<size_t> l(p.size(), 0);
    l[1] = lp[1];
    for (size_t i = 2; i < p.size(); ++i)
      l[i] = std::max(l[i - 1], lp[i]);

    return l;
  }

  /*
   * This function is to calculate the l' array for good suffix rule
   * l'[i] = largest j <= |P[i:n]| = n-i+1, such that N[j] = j
   * :param n: N array got from reverse z algorithm
   */
  std::vector<size_t> small_l_prime_array(std::vector<size_t> &n) {
    std::vector<size_t> small_lp(n.size(), 0);

    for (size_t i = 0; i < n.size(); ++i) {
      if (n[i] == i + 1) // prefix matching a suffix
        small_lp[n.size() - i - 1] = i + 1;
    }

    for (int i = n.size() - 2; i >= 0; --i) {
      if (small_lp[i] == 0) // smear them out to the left
        small_lp[i] = small_lp[i + 1];
    }

    return small_lp;
  }

  // Given a mismatch at offset i, and given L/L' and l' arrays,
  // return amount to shift as determined by good suffix rule.
  int good_suffix_mismatch(int i, std::vector<size_t> big_l_prime,
                           std::vector<size_t> small_l_prime) {
    size_t length = big_l_prime.size();
    assert(i < length);

    if (i == length - 1)
      return 1;

    i++; // i points to the leftmost matching position of P
    if (big_l_prime[i] > 0)
      return (length - big_l_prime[i]);

    return (length - small_l_prime[i]);
  }

  // calculate good suffix using border position
  // (ref:
  // https://www.geeksforgeeks.org/boyer-moore-algorithm-good-suffix-heuristic/)
  void good_suffix_rule(const std::string &pat, std::vector<size_t> &bpos,
                        std::vector<size_t> &shift) {
    size_t m = pat.size();

    size_t i = m, j = m + 1;
    bpos[i] = j;
    while (i > 0) {
      while (j <= m && pat[i - 1] != pat[j - 1]) {
        if (shift[j] == 0)
          shift[j] = j - i;
        j = bpos[j];
      }
      i--;
      j--;
      bpos[i] = j;
    }
    std::cout << "bpos list: ";
    print_arr(bpos);
    std::cout << "shift list: ";
    print_arr(shift);
  }

  void prefix_suffix_case(const std::string &pat, std::vector<size_t> &bpos,
                          std::vector<size_t> &shift) {
    size_t j = bpos.at(0);
    size_t i = 0;
    size_t m = pat.size();

    while (i <= m) {
      if (shift.at(i) == 0)
        shift.at(i) = j;
      if (i == j)
        j = bpos.at(j);
      i++;
    }
    std::cout << "final shift list: ";
    print_arr(shift);
  }

  void goodsuffix(std::string text, std::string pat) {
    size_t m = pat.size();
    size_t n = text.size();

    std::vector<size_t> shift(m + 1, 0);
    std::vector<size_t> bpos(m + 1, 0);
    good_suffix_rule(pat, bpos, shift);
    prefix_suffix_case(pat, bpos, shift);
    std::map<char, std::vector<size_t>> ex_list;
    extended_processing(pat, ex_list);

    size_t i = 0, j = m - 1;
    while (i <= n - m) {
      while (j > 0 && pat[j] != text[i + j]) {
        size_t bm_shift = get_bmshift(ex_list, text, i, j);
        size_t gs_shift = shift.at(j + 1);

        if (bm_shift > gs_shift) {
          std::cout << "choose bm_shift " << bm_shift << '\n';
          i += bm_shift;
        } else {
          std::cout << "choose gs_shift " << gs_shift << '\n';
          i += gs_shift;
        }
        j = m - 1;
      }
      j--;
      if (j < 0) {
        std::cout << "pattern matched at index " << i << " with text "
                  << text.substr(i, m) << " " << text.substr(i + m, n - i - m)
                  << '\n';
        j = m - 1;
        i += shift.at(0);
      }
    }
  }

  size_t boyer_moore(std::string p, std::string t) {
    size_t i = 0, miss = 0, match = 0;
    std::vector<size_t> match_indices;
    std::map<char, std::vector<size_t>> ex_list;
    extended_processing(p, ex_list);

    std::vector<size_t> narray = n_array(p);
    std::vector<size_t> biglprime = big_l_prime_array(p, narray);
    std::vector<size_t> bigl = big_l_array(p, biglprime);
    std::vector<size_t> smalllprime = small_l_prime_array(narray);

    size_t m = p.length(), n = t.length();
    while (i < n - m + 1) {
      size_t shift = 1;
      bool mismatch = false;
      for (int j = m - 1; j >= 0; --j) {
        if (p[j] != t[i + j]) {
          miss++;

          size_t skip_bc = get_bmshift(ex_list, t, i, j);
          size_t skip_gs = good_suffix_mismatch(j, biglprime, smalllprime);

          shift = std::max(skip_bc, skip_gs);
          // std::cout<<"\nmismatch happend, bad character shift: "<< skip_bc <<
          // " , good suffix shift: "
          //     << skip_gs << '\n';
          mismatch = true;
          // print_mismatch(p, t, i, j, shift);
          break;
        }
      }
      if (mismatch == false) {
        match++;

        match_indices.push_back(i);
        // std::cout<<"\npattern matched at index "<<i<<" with text "
        //     << t.substr(0, i) <<" " << t.substr(i, m) <<" "<< t.substr(i+m,
        //     n-i-m) << '\n';
        size_t skip_gs = smalllprime.size() - smalllprime[1];
        shift = std::max(shift, skip_gs);
      }
      i += shift;
    }
    std::cout << "all occurence indices:\n";
    print_arr(match_indices);

    return match_indices.size();
  }
};

// reads a FASTA format file line-by-line, skipping the "name" lines
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

// int main(){
//     Boyer_Moore bm = Boyer_Moore();
//     std::string T ("ACTTATTBTTCTT");
//     std::string P("CTT");
//     std::map<char, std::vector<size_t>> ex_list;
//     // bm.extended_processing(P, ex_list);
//     // bm.badchar(T, P);
//     bm.boyer_moore(P, T);
//     // std::reverse(T.begin(), T.end());
//     // std::cout<<T;

// }

int main(int argc, const char *const argv[]) {
  if (argc != 3) {
    std::cerr << "usage: " << argv[0] << " <PATTERN> <TEXT>" << std::endl;
    return EXIT_FAILURE;
  }

  const std::string P(argv[1]);
  std::string T;
  read_fasta(argv[2], T);

  // make sure pattern not bigger than text
  assert(P.length() <= T.length());

  // initialize
  Boyer_Moore bm = Boyer_Moore();

  size_t occurence = bm.boyer_moore(P, T);
  std::cout << "\ntotal number of matching occcurence is: " << occurence
            << '\n';
}
