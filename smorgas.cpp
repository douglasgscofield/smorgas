// smorgas.cpp  (c) Douglas G. Scofield, Dept. Plant Physiology, Umeå University
//
// Read pileup output and compute heterozygosity and a bunch of other things
//

// CHANGELOG
//
//
//
// TODO
// --- read pileup without indels
// --- handle indels
// --- analyze high-quality coverage
// --- raw heterozygosity
// --- model-based heterozygosity
// --- multiple pileup files
//

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
#include <stdint.h>

// #define NDEBUG  // uncomment to remove assert() code
#include <assert.h>

#include "PileupParser.h"

#include "SimpleOpt.h"

#define _SMORGAS_MAIN

#include "smorgas.h"
#include "smorgas_util.h"

using namespace std;
using namespace PileupTools;
using namespace smorgas;

static string       input_file;
static string       output_file;
static bool         opt_stdio = false;
static bool         opt_mappingquality = false;
static bool         opt_profile = false;
#ifdef _WITH_DEBUG
static int32_t      opt_debug = 1;
static int32_t      debug_progress = 100000;
static int64_t      opt_reads = -1;
static int64_t      opt_progress = 0; // 1000000;
#endif
static const char tab = '\t';
static const string sep = "\t";
//static int          map_quality_base = 33;
//static int          read_quality_base = 33;


//-------------------------------------


#ifdef _STANDALONE
int 
main(int argc, char* argv[]) {
    return main_smorgas(argc, argv);
}
#endif


//-------------------------------------


static int
usage(bool longer = false)
{
    cerr << endl;
    cerr << "Usage:   " << NAME << " [options] <in.pileup>" << endl;
    cerr << "\n\
Digest samtools mpileup output.\n\
\n\
NOTE: This command is very much a work in progress.\n\
\n";
    if (longer) cerr << "\
\n\
\n";
    cerr << "\
Options: -i FILE | --input FILE    input file name [default is stdin].  The\n\
                                   file name may also be specified on the\n\
                                   command line without this opiton.\n\
         -o FILE | --output FILE   output file name [default is stdout]\n\
         --mapping-quality         per-position mapping quality summary, to stdout\n\
         --profile                 convert to profile output for mlRho, to stdout\n\
         -? | --help               longer help\n\
\n";
#ifdef _WITH_DEBUG
    cerr << "\
         --debug INT      debug info level INT [" << opt_debug << "]\n\
         --reads INT      only process INT reads [" << opt_reads << "]\n\
         --progress INT   print reads processed mod INT [" << opt_progress << "]\n\
\n";
#endif
    cerr << endl;

    return EXIT_FAILURE;
}


//-------------------------------------


int 
smorgas::main_smorgas(int argc, char* argv[])
{

    //----------------- Command-line options

    enum { OPT_input, OPT_output, OPT_stdio,
        OPT_mappingquality,
        OPT_profile,
        OPT_opt2, OPT_opt3, OPT_opt4,
#ifdef _WITH_DEBUG
        OPT_debug, OPT_reads, OPT_progress,
#endif
        OPT_help };

    CSimpleOpt::SOption smorgas_options[] = {
        { OPT_mappingquality,  "--mapping-quality",  SO_NONE }, 
        { OPT_profile,         "--profile",          SO_NONE }, 
        { OPT_opt2,          "--opt2",         SO_NONE }, 
        { OPT_opt3,          "--opt3",         SO_NONE }, 
        { OPT_opt4,          "--opt4",         SO_NONE }, 
        { OPT_input,         "-i",                SO_REQ_SEP },
        { OPT_input,         "--input",           SO_REQ_SEP },
        { OPT_output,        "-o",                SO_REQ_SEP },
        { OPT_output,        "--output",          SO_REQ_SEP },
        { OPT_stdio,         "-",                 SO_NONE }, 
#ifdef _WITH_DEBUG
        { OPT_debug,         "--debug",           SO_REQ_SEP },
        { OPT_reads,         "--reads",           SO_REQ_SEP },
        { OPT_progress,      "--progress",        SO_REQ_SEP },
#endif
        { OPT_help,          "--help",            SO_NONE },
        { OPT_help,          "-?",                SO_NONE }, 
        SO_END_OF_OPTIONS
    };

    CSimpleOpt args(argc, argv, smorgas_options);

    while (args.Next()) {
        if (args.LastError() != SO_SUCCESS) {
            cerr << NAME << " invalid argument '" << args.OptionText() << "'" << endl;
            return usage();
        }
        if (args.OptionId() == OPT_help) {
            return usage(true);
        } else if (args.OptionId() == OPT_input) {
            input_file = args.OptionArg();
        } else if (args.OptionId() == OPT_output) {
            output_file = args.OptionArg();
        } else if (args.OptionId() == OPT_stdio) {
            opt_stdio = true;
        } else if (args.OptionId() == OPT_mappingquality) {
            opt_mappingquality = true;
        } else if (args.OptionId() == OPT_profile) {
            opt_profile = true;
#ifdef _WITH_DEBUG
        } else if (args.OptionId() == OPT_debug) {
            opt_debug = args.OptionArg() ? atoi(args.OptionArg()) : opt_debug;
        } else if (args.OptionId() == OPT_reads) {
            opt_reads = strtoll(args.OptionArg(), NULL, 10);
        } else if (args.OptionId() == OPT_progress) {
            opt_progress = args.OptionArg() ? strtoll(args.OptionArg(), NULL, 10) : opt_progress;
#endif
        } else {
            cerr << NAME << " unprocessed argument '" << args.OptionText() << "'" << endl;
            return EXIT_FAILURE;
        }
    }

    if (DEBUG(1) and ! opt_progress)
        opt_progress = debug_progress;

    if (input_file.empty()) {
        if (args.FileCount() > 1) {
            cerr << NAME << " requires at most one pileup file specified as input" << endl;
            return usage();
        } else if (args.FileCount() == 1) input_file = args.File(0);
        else input_file = "/dev/stdin";
    }

    if (output_file.empty()) {
        output_file = "/dev/stdout";
    }


    //-----------------


    PileupParser  parser(input_file);
    parser.min_base_quality = 66;
    parser.min_map_quality = 33;
    parser.debug_level = 1;

    // TODO: multiple samples
    // we will eventually handle multiple samples, for now assume whole pile is
    // described in one set of call and quality (+ mapping quality) columns
#define PRINT_UCHAR(__c__) __c__ << ":" << static_cast<uint16_t>(__c__)

    // print out mapping quality, coverage, and high-quality coverage summary per position
    if (0) {
        while (parser.read_line()) {
            parser.parse_line();
            vector<uchar_t> map_q(parser.pileup.get_map_q());
            uchar_t min_map_q = *min_element(map_q.begin(), map_q.end());
            uchar_t max_map_q = *max_element(map_q.begin(), map_q.end());
            cout << setw(8) << parser.NL << ":";
            cout << " map_q=[" << uint16_t(min_map_q) << "," << uint16_t(max_map_q) << "]";
            cout << " cov=" << parser.pileup.cov;
            size_t n_map_max = 0; uchar_t high_map_q = 58;
            for (size_t i = 0; i < map_q.size(); ++i) 
                if (map_q[i] >= high_map_q) ++n_map_max;
            cout << " cov_hiqh_q=" << n_map_max;
            cout << " frac=" << (parser.pileup.cov ? (float(n_map_max) / parser.pileup.cov) : 0);
            cout << endl;
        }
    }

    // print per-position profile for mlRho
    if (opt_profile) {
        parser.debug_level = 0;
        string current_reference = "";
        while (parser.read_line()) {
            parser.parse_line();
            if (parser.pileup.ref != current_reference) {
                // new reference
                current_reference = parser.pileup.ref;
                cout << ">" << current_reference << endl;
            }
            BaseCount bc = parser.pileup.base_count();
            cout << parser.pileup.pos;
            cout << tab << bc['A'] << tab << bc['C'] << tab << bc['G'] << tab << bc['T'] << endl;
        }
    }

    // print per-position mapping quality summary
    if (opt_mappingquality) {
        parser.debug_level = 0;
        cout << "#ref";
        cout << tab << "pos";
        cout << tab << "cov";
        cout << tab << "mapq0";
        cout << tab << "mapq60";
        cout << endl;
        while (parser.read_line()) {
            parser.parse_line();
            uchar_t mapq_a = 0; size_t mapq_a_count = 0;
            uchar_t mapq_b = 60; size_t mapq_b_count = 0;
            if (parser.pileup.cov > 0) {
                for (Pile::const_iterator citer = parser.pileup.pile.begin();
                    citer != parser.pileup.pile.end(); ++citer) {
                    uchar_t mapq = citer->map_q - parser.min_map_quality;
                    if (mapq == mapq_b) ++mapq_b_count;
                    else if (mapq == mapq_a) ++mapq_a_count;
                }
            }
            cout << parser.pileup.ref;
            cout << tab << parser.pileup.pos;
            cout << tab << parser.pileup.cov;
            cout << tab << mapq_a_count;
            cout << tab << mapq_b_count;
            cout << endl;
        }
    }

    //cout << "range base qual seen:\t" << PRINT_UCHAR(parser.min_base_quality_seen) << "\t" 
    //    << PRINT_UCHAR(parser.max_base_quality_seen) << endl;
    //cout << "range map qual seen:\t" << PRINT_UCHAR(parser.min_map_quality_seen) << "\t"
    //    << PRINT_UCHAR(parser.max_map_quality_seen) << endl;


    parser.close();

    return EXIT_SUCCESS;
}

