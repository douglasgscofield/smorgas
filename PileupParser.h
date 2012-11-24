// PileupParser.h  (c) Douglas G. Scofield, Dept. Plant Physiology, Umeå University
//
// A class that parses samtools (m)pileup
//
// Some definitions that are used by the class are provided in the namespace
// outside the class.
//
// Some terminology:
// stratum:  the information describing a single read's contribution at a position; 
//           this includes base or indel information along with quality and mapping
//           quality information
// pile:     the total collection of strata at a position; the size of the pile
//           equals the coverage (minus effects of indels, yes...?)
// pileup:   the total information contained in a single line of pileup input; 
//           includes the reference, reference base, position, coverage along with
//           everything describing the pile
// position: the location (reference + position within reference) that this pileup 
//           describes

#ifndef _PILEUPPARSER_H_
#define _PILEUPPARSER_H_

// Std C/C++ includes
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <string>
#include <memory>
#include <limits>
#include <algorithm>
#include <ctype.h>

// #define NDEBUG  // uncomment to remove assert() code
#include <assert.h>

namespace PileupTools {

//--------------------- types and type-related utilities

typedef unsigned char uchar_t;

inline uchar_t toUchar(int x) {
    return static_cast<uchar_t>(x);
}

enum readdir_t       { RD_NONE, RD_fwd, RD_rev };
enum readstructure_t { RS_NONE, RS_start, RS_end, RS_gap };
enum indel_t         { IN_NONE, IN_ins, IN_del };

//--------------------- inline functions

inline bool isBase(uchar_t c) {
    switch (std::toupper(c)) {
        case 'A': case 'C': case 'G': case 'T': case 'N':
            return true; break;
        default: 
            return false; break;
    }
}
inline bool isBaseForward(uchar_t c) {
    switch (c) {
        case 'A': case 'C': case 'G': case 'T': case 'N':
            return true; break;
        default: 
            return false; break;
    }
}
inline bool isBaseReverse(uchar_t c) {
    switch (c) {
        case 'a': case 'c': case 't': case 'g': case 'n':
            return true; break;
        default: 
            return false; break;
    }
}
inline bool isRefDirection(uchar_t c) {
    return(c == '.' or c == ',');
}
inline bool isReadBoundary(uchar_t c) {
    return(c == '^' or c == '$');
}
inline bool isIndel(uchar_t c) {
    return(c == '+' or c == '-');
}

//--------------------- utility variables and functions

inline int32_t extractNumber(const std::string& s, size_t start, size_t& end) {
    // a bit like strtol() but with a crude check for overflow
    int32_t powers_of_10[] = {     1,      10,      100,      1000,      10000, 
                              100000, 1000000, 10000000, 100000000, 1000000000};
    //size_t  max_powers_of_10 = sizeof(powers_of_10) / sizeof(int32_t);
    size_t max_pow10 = 8;  // powers from 10^0 to 10^8, missing upper range of int32_t
    int32_t ans = 0;
    int32_t sign = +1;
    if (s[start] == '+' or s[start] == '-') { sign = (s[start] == '-') ? -1 : +1; ++start; }
    if (end == 0) {
        end = start;
        for (uchar_t c = s[end]; std::isdigit(c); c = s[++end]);
    }
    if ((end - start) > max_pow10 or ! std::isdigit(s[start])) {
        std::cerr << "extractNumber: no number or out of range " << s.substr(start, end - start + 1) << std::endl;
        return 0;
    }
    for (; start < end; ++start) ans += (powers_of_10[end - start - 1] * (s[start] - '0'));
    return (sign * ans);
}

inline uchar_t lookAhead(const std::string& s, size_t lookpos) {
    return((lookpos < s.length()) ? s[lookpos] : 0);
}

inline std::string toUpper(const std::string &s) {
    // from http://stackoverflow.com/questions/735204/convert-a-string-in-c-to-upper-case
    std::string result;
    result.reserve(s.size());
    std::transform(s.begin(), s.end(), std::back_inserter(result), 
                   static_cast<int(*)(int)>(std::toupper));
    return result;
}

inline std::string toLower(const std::string &s) {
    // from http://stackoverflow.com/questions/735204/convert-a-string-in-c-to-upper-case
    std::string result;
    result.reserve(s.size());
    std::transform(s.begin(), s.end(), std::back_inserter(result), 
                   static_cast<int(*)(int)>(std::tolower));
    return result;
}


//---------------------------------------------------------------
//--------------------- Indel class


class Indel {
public:
    Indel(const int32_t sz, const std::string& sq, const size_t strat = 0, const uchar_t mq = 0);
    Indel();
    ~Indel();

    indel_t     type;
    readdir_t   dir;
    int32_t     size;
    std::string seq;
    size_t      stratum;
    uchar_t     map_q;

    void                 print(std::ostream& os = std::cerr, const std::string sep = " ") const;
    void                 print_compact(std::ostream& os = std::cerr) const;
    const std::string    seq_qualified() const;

    friend std::ostream& operator<<(std::ostream& os, const Indel& indel);
};


//---------------------------------------------------------------
//--------------------- Stratum class


class Stratum {
public:
    Stratum();
    ~Stratum();

    uchar_t                 base;
    uchar_t                 base_q;
    uchar_t                 map_q;
    readdir_t               dir;
    readstructure_t         read_str;  // TODO: infer read mapping quality from read structure
    uchar_t                 read_map_q;
    Indel*                  indel;  // if there's an indel, this points to an element of pileup.indels[]
};
typedef std::vector<Stratum> Pile;


//---------------------------------------------------------------
//--------------------- Pileup class


class Pileup {  // one entry for each unique position
    // TODO: methods in PileupParser will walk this for analysis... better solution
    // might be to add a few methods here to assist?  or make a derived class that
    // adds those methods?
public:
    Pileup(uchar_t min_base_qual = '\0');
    ~Pileup();

    std::string             ref;
    size_t                  pos;
    uchar_t                 refbase;  // TODO: can reference base be more than one character?
    int32_t                 cov;
    // TODO: multiple samples
    Pile                    pile;  // the pile has 1+ strata TODO: is 0 ever true?
    std::vector<Indel>      indels;  // less space to keep them here and not in Stratum

    enum parsestate_t { PS_NONE=0x0, PS_lite=0x1, PS_pile=0x2, PS_all=0x3 };
    parsestate_t            parse_state;

    // TODO: allow different qualities for different samples
    uchar_t                 min_set_base_quality;
    uchar_t                 min_set_map_quality;

    bool                    set_min_base_quality(uchar_t min_base_q);
    bool                    set_min_map_quality(uchar_t min_map_q);
    std::vector<uchar_t>    get_map_q(const size_t start = 0, size_t end = std::numeric_limits<std::size_t>::max());

    void        print(std::ostream& os = std::cerr) const;
    void        print_pile(std::ostream& os = std::cerr, const size_t start = 0, size_t end = 0, 
                           const std::string sep = "\t") const;
    void        print_pile_stack(std::ostream& os = std::cerr,
                           const size_t start = 0, size_t end = 0, const bool include_indels = false, 
                           const std::string end_stack = "\n") const;
    void        debug_print();
    void        debug_print_pile(const size_t start = 0, size_t end = 0, bool stack = false);
                               
    friend std::ostream& operator<<(std::ostream& os, const Pileup& pileup);

};


//---------------------------------------------------------------
//--------------------- PileupParser class


class PileupParser {

public:
    static const std::string name()    { return "PileupParser"; }
    static const std::string version() { return "0.0.1-dev"; }
    static const std::string author()  { return "Douglas G. Scofield"; }
    static const std::string contact() { return "douglasgscofield@gmail.com"; }

public:
    std::ifstream       stream;  // stream we're reading from
    std::string         filename; // filename, if one was given
    const char          FS;      // input field separator
    const char          RS;      // input line separator
    size_t              NL;      // line number within pileup file
    int                 NF;      // number of fields in current line
    uchar_t             min_base_quality;
    uchar_t             min_map_quality;

public:

    std::string              line;  // mpileup line we're currently working on

    std::vector<std::string> fields;  // fields of mpileup line
    enum { F_ref=0, F_pos, F_refbase, F_cov, F_base_call, F_base_q, F_map_q, F_END };

    std::vector<std::string> references;  // reference sequences named in the pileup

    Pileup                   pileup;

    uchar_t     min_base_quality_seen;
    uchar_t     max_base_quality_seen;
    uchar_t     min_map_quality_seen;
    uchar_t     max_map_quality_seen;

    // ctor, dtor
    PileupParser(const std::string& fname);
    PileupParser();
    ~PileupParser();

    void        open(const std::string& fname);
    void        close();
    int         read_line();
    void        parse_line();
    void        parse_line_lite();
    void        parse_pile();

    void        print(std::ostream& os = std::cerr) const;
    void        print_lite(std::ostream& os = std::cerr, const std::string sep = "\t") const;
    void        scan(size_t n_lines);
    friend std::ostream& operator<<(std::ostream& os, const PileupParser& parser);

    void        update_qualities_seen();

    int         debug_level;
    inline bool debug(int level) { return(debug_level >= level); }


};  // class PileupParser


} // namespace PileupTools


#endif // _PILEUPPARSER_H_
