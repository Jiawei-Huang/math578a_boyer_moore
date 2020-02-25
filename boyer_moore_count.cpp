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

  //     1. bad character heuristic

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

    // std::cout << "character : index\n";
    // print_map(ex_list);
  }

  /**
   * @brief             calculate the shift space for bad character rule
   *
   * @param ex_list     extended list containing all occurence of each character
   *                    in pattern.
   * @param text        text
   * @param i           start index of the text
   * @param j           current index of the pattern
   * @return size_t     shift spaces
   */
  size_t get_bmshift(std::map<char, std::vector<size_t>> &ex_list,
                     const char &text, const size_t &i, const size_t &j) {
    std::map<char, std::vector<size_t>>::iterator it = ex_list.find(text);
    size_t shift = 0;
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
    // size_t m = tab[(int)text].size() - 1;
    // size_t i = 0;
    // while (i <= m && j < tab[(int)text][i]) {
    //  i++;
    //}
    // if (i == m + 1)
    //  return shift;
    // shift = j - tab[(int)text][i];
    return shift;
  }

  /**
   * @brief bad character boyer-moore algorithm
   *
   * @param text      text
   * @param pat       pattern
   */
  void badchar(const std::string &text, const std::string &pat) {
    size_t n = text.size();
    size_t m = pat.size();
    std::map<char, size_t> alphabet;

    processing(pat, alphabet);
    std::map<char, std::vector<size_t>> ex_list;
    extended_processing(pat, ex_list);

    size_t i = 0;
    int j = m - 1;
    while (i <= n - m) {
      while (pat[j] != text[i + j]) {
        std::cout << "mismatch happend \n";
        size_t shift = get_bmshift(ex_list, text[i + j], i, j);

        // int shift = std::max(j - alphabet.at(text[i+j]), 1);

        print_mismatch(pat, text, i, j, shift);
        i += shift;
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

  /**
   * @brief     Get the alphabet map object
   *
   * @param alphabet    alphabet we define, like "ACGT"
   * @return std::map<char, size_t>     return a map from alphabet characters to
   * integers
   */
  std::map<char, size_t> get_alphabet_map(const std::string alphabet = "ACGT") {
    std::map<char, size_t> amap;
    for (size_t i = 0; i < alphabet.size(); ++i) {
      amap[alphabet[i]] = i;
    }

    return amap;
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
    // table.resize(256);

    // size_t i = p.size();
    // while (i--) {
    //  table[(int)p[i]].push_back(i);
    //}

    for (size_t i = 0; i < p.size(); ++i) {
      char c = p[i];
      table.push_back(index);
      index[dna_encoding[c]] = i + 1;
    }

    return table;
  }

  //     2. goof suffix rule
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
    // then for any i < j, l'[i] = l'[j], here
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
    // assert(i < length);

    if (i == length - 1)
      return 1;

    i++; // i points to the leftmost matching position of P
    if (big_l_prime[i] > 0)
      return (length - big_l_prime[i]);

    return (length - small_l_prime[i]);
  }

  /**
   * @brief         calculate good suffix using border position
   *                (ref:https://www.geeksforgeeks.org/boyer-moore-algorithm-good-suffix-heuristic/)
   * @param pat     pattern
   * @param bpos        border position array, length is m+1
   * @param shift       shift array, length is m+1
   */
  void good_suffix_rule(const std::string &pat, std::vector<size_t> &bpos,
                        std::vector<size_t> &shift) {
    size_t m = pat.size();

    // start at m, backwards
    size_t i = m, j = m + 1;
    bpos[i] = j;
    while (i > 0) {
      while (j <= m && pat[i - 1] != pat[j - 1]) {
        // mismatch at i-1. shift at i
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

  /**
   * @brief     case 2, see refrence above
   *
   * @param pat     pattern
   * @param bpos    border position array get from good_suffix_rule function
   * @param shift   shift array get from good_suffix_rule function
   */
  void prefix_suffix_case(const std::string &pat,
                          const std::vector<size_t> &bpos,
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

  /**
   * @brief       boyer moore based on border position array
   *
   * @param text      text
   * @param pat       pattern`
   */
  void goodsuffix(const std::string &text, const std::string &pat) {
    size_t m = pat.size();
    size_t n = text.size();

    std::vector<size_t> shift(m + 1, 0);
    std::vector<size_t> bpos(m + 1, 0);
    good_suffix_rule(pat, bpos, shift);
    prefix_suffix_case(pat, bpos, shift);
    std::map<char, std::vector<size_t>> ex_list;
    extended_processing(pat, ex_list);

    int i = 0, j = m - 1;
    while (i <= n - m) {
      while (j > 0 && pat[j] != text[i + j]) {
        size_t bm_shift = get_bmshift(ex_list, text[i + j], i, j);
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

  /**
   * @brief       boyer moore based on Gusfield
   *
   * @param p         pattern
   * @param t         text
   * @return size_t       number of matching occurrence
   */
  void boyer_moore(const std::string &p, const std::string &t) {
    size_t i = 0, match = 0, comparision = 0;
    // std::vector<size_t> match_indices;
    // std::map<char, std::vector<size_t>> ex_list;
    // extended_processing(p, ex_list);
    std::vector<std::vector<size_t>> tab = get_bad_char_table(p);

    std::vector<size_t> narray = n_array(p);
    std::vector<size_t> biglprime = big_l_prime_array(p, narray);
    std::vector<size_t> bigl = big_l_array(p, biglprime);
    std::vector<size_t> smalllprime = small_l_prime_array(narray);

    size_t m = p.length(), n = t.length();
    auto start = std::chrono::steady_clock::now();
    while (i < n - m + 1) {
      size_t shift = 1;

      int j = m - 1;
      while (j >= 0 && p[j] == t[i + j]) {
        j--;
      }
      comparision += m - j;
      if (j < 0) {
        match++;
        comparision--;
        shift = m - smalllprime[1];
        //  size_t skip_gs = smalllprime.size() - smalllprime[1];
        //  shift = std::max(shift, skip_gs);
      } else {
        size_t skip_bc = bad_character_mismatch(tab, t[i + j], j);
        size_t skip_gs = good_suffix_mismatch(j, biglprime, smalllprime);
        shift = std::max(skip_bc, skip_gs);
      }
      i += shift;
      //   while (j--) {
      //     if (p[j] != t[i + j]) {

      //       // size_t skip_bc = get_bmshift(ex_list, t[i + j], i, j);
      //       size_t skip_bc = bad_character_mismatch(tab, t[i + j], j);
      //       size_t skip_gs = good_suffix_mismatch(j, biglprime, smalllprime);
      //       shift = std::max(skip_bc, skip_gs);
      // // print mismatch position
      // std::cout << "\nmismatch happend, bad character shift: " << skip_bc
      //           << " , good suffix shift: " << skip_gs << '\n';
      // print_mismatch(p, t, i, j, shift);
      //       mismatch = true;
      //       comparision += m - j;
      //       break;
      //     }
      //   }
      //   if (mismatch == false) {
      //     match++;
      // print matching position
      // std::cout << "\npattern matched at index " << i << " with text "
      //           << t.substr(0, i) << " " << t.substr(i, m) << " "
      //           << t.substr(i + m, n - i - m) << '\n';
      //     size_t skip_gs = smalllprime.size() - smalllprime[1];
      //     shift = std::max(shift, skip_gs);
      //     comparision += m;
      //   }
      //   i += shift;
    }
    auto end = std::chrono::steady_clock::now();
    auto diff = end - start;

    // std::cout << "all occurence indices:\n";
    // print_arr(match_indices);

    std::cout << "\ntotal number of matching: " << match << '\n';
    std::cout << "\ntotal number of comparisions: " << comparision << '\n';
    std::cout << "\ntime used for comparision: ";
    std::cout << std::chrono::duration_cast<std::chrono::seconds>(diff).count()
              << " s" << std::endl;
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

// This function just removes the sequence (e.g. chromosome) names and
// newlines from a string loaded from a FASTA file. It is not
// optimized for speed, but it should be pretty clear.
static void remove_names_newlines(std::string &T) {
  bool outside_name = true;
  size_t j = 0;
  const size_t n = T.size();
  for (size_t i = 0; i < n; ++i) {
    const char c = T[i];
    if (outside_name) {
      if (c == '>')
        outside_name = false;
      else if (c != '\n') {
        T[j++] = c;
      }
    } else
      outside_name = (c == '\n');
  }
  // resize but keep capacity
  T.resize(j);
}

// The "get_filesize" below uses some pretty specific C++ code for
// "streams" and if you don't understand it, that's fine.
static size_t get_filesize(const std::string &filename) {
  std::ifstream in(filename);
  if (!in)
    throw std::runtime_error("problem with file: " + filename);
  const std::streampos begin_pos = in.tellg();
  in.seekg(0, std::ios_base::end);
  return in.tellg() - begin_pos;
}

static void read_fasta_as_one_sequence(const std::string &fasta_filename,
                                       std::string &T) {
  T.clear(); // start with empty string

  const size_t filesize = get_filesize(fasta_filename);

  // using "C" functions to read in the input because on my mac there
  // is a problem with reading large files using a single "read"
  // function call when using C++ streams...
  FILE *in = fopen(fasta_filename.c_str(), "rb");
  if (!in)
    throw std::runtime_error("problem with file: " + fasta_filename);

  T.resize(filesize); // change *size*, not capacity of T here

  if (fread((char *)&T[0], 1, filesize, in) != filesize)
    throw std::runtime_error("problem with file: " + fasta_filename);

  if (fclose(in) != 0)
    throw std::runtime_error("problem with file: " + fasta_filename);

  // remove the sequence names from the FASTA format string, along
  // with the newline characters, what remains should be just DNA
  // bases (maybe with a few random IUPAC degenerate nucleotides)
  remove_names_newlines(T);
}

int main(int argc, const char *const argv[]) {
  if (argc != 3) {
    std::cerr << "usage: " << argv[0] << " <PATTERN> <TEXT>" << std::endl;
    return EXIT_FAILURE;
  }

  const std::string P(argv[1]);
  std::string T;
  read_fasta_as_one_sequence(argv[2], T);

  // make sure pattern not bigger than text
  assert(P.length() <= T.length());

  // initialize
  Boyer_Moore bm = Boyer_Moore();

  bm.boyer_moore(P, T);
}
