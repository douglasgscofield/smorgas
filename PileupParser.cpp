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

#include "PileupParser.h"

namespace PileupTools {


//--------------------------------------------------------
//--------------------------------- class PileupParser


//----------------- ctor and dtor

PileupParser::PileupParser(const std::string& fname)
    : filename(fname), FS('\t'), RS('\n'), NL(0), NF(0), 
      minimum_base_quality(0), minimum_map_quality(0), line(""), fields(7),
      minimum_base_quality_character_seen(0xff), maximum_base_quality_character_seen(0x00),
      minimum_map_quality_character_seen(0xff), maximum_map_quality_character_seen(0x00),
      debug_level(0)
{
    open(filename);
}

PileupParser::PileupParser()
    : filename(""), FS('\t'), RS('\n'), NL(0), NF(0), 
      minimum_base_quality(0), minimum_map_quality(0), line(""), fields(7),
      minimum_base_quality_character_seen(0xff), maximum_base_quality_character_seen(0x00),
      minimum_map_quality_character_seen(0xff), maximum_map_quality_character_seen(0x00),
      debug_level(0)
{
}

PileupParser::~PileupParser()
{
}

//----------------- file handling and raw line/field reading

void 
PileupParser::open(const std::string& fname)
{
    filename = fname;
    stream.open(filename.c_str());
    NL = 0;
    line = "";
}

void 
PileupParser::close()
{
    stream.close();
}

int 
PileupParser::read_line()
{
    NF = 0;
    if (getline(stream, line, RS)) {
    	++NL;
        if (debug(2)) std::cerr << "line " << NL << " :" << line << ":" << std::endl;
        // process line fields
        size_t pos = 0, f;
        for (f = 0; f < F_END; ++f) {
            size_t t = line.find(FS, pos);
            if (t != std::string::npos) {
                fields[f] = line.substr(pos, t - pos);
                if (debug(3)) std::cerr << "field " << f << " :" << fields[f] << ":" << std::endl;
                pos = t + 1;  // skip the FS
            } else {  // no FS, so last field is the remainder of the line
                fields[f] = line.substr(pos);
                if (debug(3)) std::cerr << "last field " << f << " :" << fields[f] << ":" << std::endl;
                break;
            }
        }
        NF = f + 1;
    }
    return(NF);
}

//----------------- parse the current line/fields

void
PileupParser::parse_line()
{
    const char* const thisfunc = "parse_line";
    if (line == "") { std::cerr << thisfunc << ": no line to parse" << std::endl; return; }
    parse_line_lite();
    parse_pile();
    // here, parse_state will be (PS_lite | PS_pile) == PS_all
}

void
PileupParser::parse_line_lite()
{
    const char* const thisfunc = "parse_line_lite";
    if (line == "") { std::cerr << thisfunc << ": no line to parse" << std::endl; return; }
    if (references.size() == 0 or fields[F_ref] != references.back()) { // assumes input sorted
        references.push_back(fields[F_ref]);
    }
    pileup.ref = references.back();
    pileup.pos = atol(fields[F_pos].c_str());
    pileup.refbase = fields[F_refbase][0];
    pileup.cov = atol(fields[F_cov].c_str());
    pileup.parse_state = static_cast<Pileup::parsestate_t>(pileup.parse_state | Pileup::PS_lite);
    // do not do parse_pile() here
}

void
PileupParser::parse_pile()
{
    const char* const thisfunc = "parse_pile";

    Pile& pile = pileup.pile;

    pile.resize(pileup.cov);

    size_t stratum = 0; // position within pile (in terms of strata)
    size_t i = 0; // position within base string, contains other info

    while (i < fields[F_base_call].length()) {

        if (stratum == pile.size()) {
            std::cerr << "needed to resize vector from " << pileup.cov << std::endl;
            pile.resize(pileup.cov + int(pileup.cov / 2));
        }

        // Each base call entry can be one of several types.  
        // TODO: does this cover it?  Can we have IUPAC or length > 1?
        //
        // B = [ACGTN]
        // b = [acgtn]
        // # = numeric character
        // q = quality character
        //
        // ^q[.,]    : read begins here, mapping quality q, forward/reverse orientation, matches ref
        // ^qB       : read begins here, mapping quality q, forward orientation, does not match ref
        // ^qb       : read begins here, mapping quality q, reverse orientation, does not match ref
        // [.,]      : read in forward (.) or reverse (,) orientation, matches ref
        // [.,]$     : read ends here, forward/reverse orientation, matches ref
        // .[+-]#+B+ : indel length #+, contents B+, forward orientation
        // ,[+-]#+b+ : indel length #+, contents b+, reverse orientation
        // B         : read in forward orientation, does not match ref
        // b         : read in reverse orientation, does not match ref
        // B$        : read ends here, forward orientation, does not match ref
        // b$        : read ends here, reverse orientation, does not match ref
        // *         : position is a continuation of a deletion in the read at this stratum
        //

        uchar_t c0 = fields[F_base_call][i];

        if (c0 == '^') {  // if read start, eat it and move to next character

            pile[stratum].read_str = RS_start;
            pile[stratum].read_map_q = fields[F_base_call][i + 1];
            i += 2;
            c0 = fields[F_base_call][i];

        }

        switch (c0) { // read direction (. or ,) or base, optionally followed by $, or *
            case '.': 
                pile[stratum].dir = RD_fwd;
                pile[stratum].base = pileup.refbase;
                break;
            case ',': 
                pile[stratum].dir = RD_rev;
                pile[stratum].base = pileup.refbase;
                break;
            case 'A': case 'C': case 'G': case 'T': case 'N':
                pile[stratum].dir = RD_fwd;
                pile[stratum].base = c0;
                break;
            case 'a': case 'c': case 'g': case 't': case 'n':
                pile[stratum].dir = RD_rev;
                pile[stratum].base = toupper(c0);
                break;
            case '*':
                pile[stratum].base = '*';
                pile[stratum].read_str = RS_gap;
                break;
            default:
                std::cerr << thisfunc << ": line " << NL << " stratum " << stratum
                    << " unknown base call character: " << c0 << std::endl;
                break;
        }

        uchar_t c1 = lookAhead(fields[F_base_call], (i + 1));  // returns next char or 0 if end of string

        if (isIndel(c1)) {   // [+-]#+[Bb]+

            // we have already seen the leading . or ,
            // eat [+-]#+ for indel size, then use abs(indel size) to eat the sequence

            size_t j = i + 1, k = 0;  // start at '+' or '-', k will be first non-numeric char
            int32_t indel_size = extractNumber(fields[F_base_call], j, k);
            Indel indel(indel_size, fields[F_base_call].substr(j, abs(indel_size)), stratum);
            pileup.indels.push_back(indel);
            pile[stratum].indel = &pileup.indels.back();  // pointer to its new spot in indels stack
            i = k - 1;

        } else if (c1 == '$') {

            pile[stratum].read_str = RS_end;
            i += 1;

        }

        // after all that mess, the base and mapping quality columns are easy
        pile[stratum].base_q = fields[F_base_q][stratum];
        if (fields[F_map_q] != "") {
            pile[stratum].map_q = fields[F_map_q][stratum];
            if (pile[stratum].indel)
                pile[stratum].indel->map_q = pile[stratum].map_q;
        }

        if (pile[stratum].read_map_q and pile[stratum].map_q
            and pile[stratum].map_q != pile[stratum].read_map_q)
            std::cerr << "line " << NL << " stratum " << stratum << " read_map_q "
                << pile[stratum].read_map_q << " != map_q " << pile[stratum].map_q 
                << std::endl;

        update_qualities_seen(pile[stratum].base_q, pile[stratum].map_q);

        ++stratum;
        ++i;
    }

    if (pile.size() != stratum) {
        std::cerr << "at end of parse_line, pile.size() " << pile.size() << " != stratum "
            << stratum << std::endl;
        if (pile.size() > stratum)
            std::cerr << "MAJOR PROBLEM: shrinking the size of the pile!!!!!!" << std::endl;
        pile.resize(stratum);
    }

    pileup.parse_state = static_cast<Pileup::parsestate_t>(pileup.parse_state | Pileup::PS_pile);

    if (debug(1)) 
        std::cerr << thisfunc << ": line " << NL << " has " << stratum << " strata" << std::endl;
    //assert(pile.size() == stratum);
}


//----------------- other stuff

void
PileupParser::scan(size_t n_lines)
{
    // scan the input pileup and do a quick summary of what is seen
    // default to say 1000 lines, on option do more up to whole file
}


//--------------------------------------------------------
//--------------------------------- class Pileup


Pileup::Pileup(uchar_t min_base_qual, uchar_t min_map_qual) 
    : ref(""), pos(0), refbase('\0'), cov(-1), parse_state(PS_NONE), 
      minimum_base_quality(min_base_qual), minimum_map_quality(min_map_qual)
{
}


Pileup::Pileup() 
    : ref(""), pos(0), refbase('\0'), cov(-1), parse_state(PS_NONE),
      minimum_base_quality('\0'), minimum_map_quality('\0')
{
}


Pileup::~Pileup()
{
    // when out of scope, we get automatic delete of unique_ptr<Pile> pile
}
   

void 
Pileup::print(std::ostream& os) const
{
    os << "PS=" << parse_state;
    os << " ref=" << ref;
    os << " pos=" << pos;
    os << " refbase=" << refbase;
    os << " cov=" << cov;
    os << " min-base-q=" << minimum_base_quality;
    os << " min-map-q=" << minimum_map_quality;
    if (pile.size() == 0) {
        os << " NO_PILE "; 
    } else {
        os << " base="; 
        for (size_t i = 0; i < pile.size(); ++i) {
            os << (pile[i].dir == RD_fwd ? pile[i].base : tolower(pile[i].base));
            if (pile[i].indel) indel->print_compact(os);
        }
        os << " base_q="; for (size_t i = 0; i < pile.size(); ++i) os << pile[i].base_q;
        os << " map_q="; for (size_t i = 0; i < pile.size(); ++i) os << pile[i].map_q;
    }
}


std::ostream& operator<<(std::ostream&os, const Pileup& pileup)
{
    os << "pileup:";
    pileup.print(os);
    return(os);
}


bool 
Pileup::apply_minimum_base_quality(uchar_t min_base_q)
{
    if (min_base_q == 0)
        std::cerr << "pileup::apply_minimum_base_quality: argument is 0, so a no-op" << std::endl;
    bool good_quals = true;
    if (minimum_base_quality > 0)
        std::cerr << "Pileup::apply_minimum_base_quality: minimum base quality already nonzero" << std::endl;
    minimum_base_quality = min_base_q;
    for (size_t i = 0; i < pile.size(); ++i) {
        if (pile[i].base_q < min_base_q)
            good_quals = false;
        pile[i].base_q -= min_base_q;
    }
    return good_quals;
}


bool 
Pileup::apply_minimum_map_quality(uchar_t min_map_q)
{
    if (min_map_q == 0)
        std::cerr << "pileup::apply_minimum_base_quality: argument is 0, so a no-op" << std::endl;
    bool good_quals = true;
    if (minimum_map_quality > 0)
        std::cerr << "pileup::apply_minimum_base_quality: minimum map quality already nonzero" << std::endl;
    minimum_map_quality = min_map_q;
    for (size_t i = 0; i < pile.size(); ++i) {
        if (pile[i].map_q < min_map_q or pile[i].read_map_q < min_map_q)
            good_quals = false;
        pile[i].map_q -= min_map_q;
        pile[i].read_map_q -= min_map_q;
    }
    return good_quals;
}


//--------------------------------------------------------
//--------------------------------- class Stratum


Stratum::Stratum() 
    : base('\0'), base_q('\0'), map_q('\0'), dir(RD_NONE), read_str(RS_NONE),
      read_map_q('\0'), indel(0)
{
}

Stratum::~Stratum()
{
}


//--------------------------------------------------------
//--------------------------------- class Indel


Indel::Indel(const int32_t sz, const std::string& sq, const size_t strat = 0, const uchar_t mq = 0)
    : type(sz > 0 ? IN_ins : IN_del),
      dir(isBaseForward(sq[0]) ? RD_fwd : RD_rev),
      size(sz),
      seq(dir == RD_fwd ? sq : toUpper(sq)),
      stratum(strat),
      map_q(mq)
{
    // all the business is taken care of with initializers
    // sz: signed, with + = insertion, - = deletion
    // sq: case-sensitive, with upper = forward-direction, lower = reverse-direction
    // strat: the stratum to which this indel belongs
    // TODO: maybe we want to add base qualities and mapping quality of the read?
}


Indel::Indel()
    : type(IN_NONE), dir(RD_NONE), size(0), seq(""), stratum(0), map_q(0)
{
}


Indel::~Indel()
{
}


void
Indel::print(std::ostream& os) const
{
    os << "(@" << stratum << " ";
    if (type == IN_NONE)
        os << "0)";
    else
        os << (dir == RD_fwd ? "." : ",") << " " << size << " " << seq_qualified() << ")";
}


void
Indel::print_compact(std::ostream& os) const
{
    if (type == IN_NONE)
        os << "(0)";
    else
        os << "(" << seq_qualified() << ")";
}


const std::string&
Indel::seq_qualified() const
{
    std::ostreamstream qseq;
    qseq << (type > IN_ins ? "+" : "-") << (dir == RD_rev ? toLower(seq) : seq);
    return(qseq.str())
}


std::ostream& operator<<(std::ostream&os, const Indel indel)
{
    os << "indel";
    indel.print(os);
    return(os);
}


} // namespace PileupTools
