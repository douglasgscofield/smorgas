// smorgas.cpp  (c) Douglas G. Scofield, Dept. Plant Physiology, Umeå University
//
// Read pileup output and compute heterozygosity and a bunch of other things
//

// CHANGELOG
//
//
//
// TODO
// -x- read pileup without indels
// -x- handle indels
// --- base qualities (base 64) can be < 64 in vicinity of N, and ! (33) when the reference is N
/* 
MA_1    1098    A       15      ,,,,,,,,,,,,,,. \\\\\\Z\\\\\Y\\
MA_1    1099    T       15      ,,,,,,,,,,,,,,. aaaaaa`aaaa^Taa
MA_1    1100    T       15      ,,,,,,,,,,,,,,. aaaa_aYaaaaTYaa
MA_1    1101    G       15      ,,,,,,,,,,,,,,. \\\\\\\\\\\\\\\
MA_1    1102    G       15      ,,,,,,,,,,,,,,. ZZZZZZZZZZZZZZZ
MA_1    1103    G       15      ,,,,,,,,,,,,,,. ZZZZZZZZZZZZZZZ
MA_1    1104    G       15      ,,,,,,,,,,,,,,. YYYYY\Y[YYYYYY\
MA_1    1105    T       15      ,,,,,,,,,,,,,,. ?????>?>???????
MA_1    1106    T       15      ,,,,,,,,$,,,,,,.        <<<<<;<;<<<<<<<
MA_1    1107    T       14      ,-2nn,-2nn,-2nn,-2nn,-2nn,$,-2nn,-2nn,-2nn,-2nn,-2nn,-2nn,-2nn. :::::9::::::::
MA_1    1108    N       13      ************G   !!!!!!!!!!!!9
MA_1    1109    N       13      ************G   !!!!!!!!!!!!8
MA_1    1110    N       13      ggggggggggggA$  !!!!!!!!!!!!7
MA_1    1111    N       12      gggggggggggg    !!!!!!!!!!!!
MA_1    1112    N       12      aaaaaaaaaaaa    !!!!!!!!!!!!
MA_1    1113    N       12      aaaaaaaaaaaa    !!!!!!!!!!!!
MA_1    1114    N       12      c$c$c$c$c$c$c$c$c$c$c$c$        !!!!!!!!!!!!
MA_1    1904    N       1       ^Fa     5
MA_1    1905    N       1       t       5
MA_1    1906    N       1       g       6
MA_1    1907    N       4       ^3t^1t^3tt      7776
MA_1    1908    N       5       ^3ggggg 87777
MA_1    1909    N       6       aaaa^3aa        888898
MA_1    1910    N       7       cccccc^FC       999999;
MA_1    1911    N       7       ggggggG ;;;;;;;
MA_1    1912    G       14      ,^>.^6.,^0.,,,^>.^6.^6.^>.,.    >>>>>>>>>>>>>>
MA_1    1913    T       15      ,..,.,,,....^8.,.       HHHHHHHHHHHHBHH
*/
// --- analyze high-quality coverage
// --- raw heterozygosity
// --- model-based heterozygosity
// --- multiple BAMs in pileup input; handle with stack that parallels the strata
//     and indicates which BAM in sequence they came from... or maybe keep track
//     only of intervals that give us slices from which to take raw and parsed
//     data
// --- check for unsorted input
// --- solve read mapping quality from reads if separate column not available;
//     handle with stack that parallels the strata and tracks read quality for
//     each stratum as reads are introduced and removed
//

#include "PileupParser.h"

namespace PileupTools {


//--------------------------------------------------------
//--------------------------------- class PileupParser


//----------------- ctor and dtor

PileupParser::PileupParser(const std::string& fname)
    : filename(fname), FS('\t'), RS('\n'), NL(0), NF(0), 
      min_base_quality(0), min_map_quality(33), line(""), fields(7),
      min_base_quality_seen(0xff), max_base_quality_seen(0x00),
      min_map_quality_seen(0xff), max_map_quality_seen(0x00),
      debug_level(0)
{
    open(filename);
}

PileupParser::PileupParser()
    : filename(""), FS('\t'), RS('\n'), NL(0), NF(0), 
      min_base_quality(0), min_map_quality(33), line(""), fields(7),
      min_base_quality_seen(0xff), max_base_quality_seen(0x00),
      min_map_quality_seen(0xff), max_map_quality_seen(0x00),
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
        for (f = 0; f < F_END; ++f) {  // TODO: handle multiple BAMs here
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

// parse_line_lite() should do a cursory parse of the fields and load access to them
// into the pileup object.  If there are multiple BAMs in the pileup, then here is
// where the columns should be joined and the masks that allow methods to extract info
// for each BAM are generated; the final format here for base call, base quality and
// mapping quality (if supplied), though, is a single string for each.

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
    pileup.raw_base_call = &fields[F_base_call];
    pileup.raw_base_quality = &fields[F_base_q];
    pileup.raw_map_quality = &fields[F_map_q];
    pileup.parse_state = Pileup::PS_lite;
    // do not do parse_pile() here
}

void
PileupParser::parse_pile()
{
    const char* const thisfunc = "parse_pile";

    pileup.reset_pile();  // does not affect anything done by parse_line_lite()

    Pile& pile = pileup.pile;

    pile.resize(pileup.cov);

    size_t stratum = 0; // position within pile (in terms of strata)
    size_t i = 0; // position within base string, contains other info

    while (i < pileup.raw_base_call->length()) {

        if (pileup.cov == 0) break;  // the base call column may still hold '*' so don't even go there

        if (stratum == pile.size()) {
            std::cerr << "NL=" << NL << " i=" << i <<" stratum=" << stratum 
                << " RESIZING pile from " << pile.size();
            pile.resize(pile.size() + 1);
            //pile.resize(pileup.cov + int(pileup.cov / 2));
            std::cerr << " to " << pile.size() << std::endl;
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

        uchar_t c0 = (*pileup.raw_base_call)[i];

        if (c0 == '^') {  // if read start, eat it and move to next character

            pile[stratum].read_str = RS_start;
            pile[stratum].read_map_q = (*pileup.raw_base_call)[i + 1];
            i += 2;
            c0 = (*pileup.raw_base_call)[i];

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

        uchar_t c1 = lookAhead((*pileup.raw_base_call), (i + 1));  // returns next char or 0 if end of string

        if (isIndel(c1)) {   // [+-]#+[Bb]+

            // we have already seen the leading . or ,
            // eat [+-]#+ for indel size, then use abs(indel size) to eat the sequence

            size_t j = i + 1, k = 0;  // start at '+' or '-', k will be first non-numeric char
            int32_t indel_size = extractNumber((*pileup.raw_base_call), j, k);
            Indel indel(indel_size, pileup.raw_base_call->substr(k, abs(indel_size)), stratum);
            pileup.indels.push_back(indel);
            pile[stratum].indel = &pileup.indels.back();  // new spot in indels stack
            i = k + abs(indel_size) - 1;  // i points to last char of indel sequence

        } else if (c1 == '$') {

            pile[stratum].read_str = RS_end;
            i += 1;

        }

        // after all that mess, the base and mapping quality columns are easy
        if (stratum < pileup.raw_base_quality->size()) 
            pile[stratum].base_q = (*pileup.raw_base_quality)[stratum];
        else
            std::cerr << "NL=" << NL << " stratum=" << stratum 
                << " exceeds length of base_q" << std::endl;

        if ((*pileup.raw_map_quality) != "") {
            if (stratum < pileup.raw_map_quality->size()) 
                pile[stratum].map_q = (*pileup.raw_map_quality)[stratum];
            else
                std::cerr << "NL=" << NL << " stratum=" << stratum 
                    << " exceeds length of map_q" << std::endl;
            if (pile[stratum].indel)
                pile[stratum].indel->map_q = pile[stratum].map_q;
        }

        if (pile[stratum].read_map_q and pile[stratum].map_q
            and pile[stratum].map_q != pile[stratum].read_map_q) {
            std::cerr << "NL=" << NL << " stratum=" << stratum 
                << " read_map_q != map_q: " << pile[stratum].read_map_q 
                << " vs " << pile[stratum].map_q << std::endl;
        }

        ++stratum;
        ++i;
    }

    if (debug(2)) 
        std::cerr << thisfunc << ": line " << NL << " has " << stratum << " strata" << std::endl;

    update_qualities_seen();
    //if (min_base_quality) 
    //    pileup.set_min_base_quality(min_base_quality);
    //if (min_map_quality) 
    //    pileup.set_min_map_quality(min_map_quality);

    if (pile.size() != stratum) {
        std::cerr << "at end of parse_line, pile.size() " << pile.size() << " != stratum "
            << stratum << std::endl;
        if (pile.size() > stratum)
            std::cerr << "MAJOR PROBLEM: shrinking the size of the pile!!!!!!" << std::endl;
        pile.resize(stratum);
    }

    pileup.parse_state = static_cast<Pileup::parsestate_t>(pileup.parse_state | Pileup::PS_pile);

    //assert(pile.size() == stratum);
}


//----------------- other stuff


void
PileupParser::scan(size_t n_lines)
{
    // scan the input pileup and do a quick summary of what is seen
    // default to say 1000 lines, on option do more up to whole file
}


void 
PileupParser::update_qualities_seen()
{
    uchar_t prev_min_bq = min_base_quality_seen, prev_max_bq = max_base_quality_seen, 
            prev_min_mq = min_map_quality_seen, prev_max_mq = max_map_quality_seen;
    if (pileup.cov == 0) return;
    for (size_t i = 0; i < pileup.pile.size(); ++i) {
        if (pileup.pile[i].base_q < min_base_quality_seen)
            min_base_quality_seen = pileup.pile[i].base_q;
        else if (pileup.pile[i].base_q > max_base_quality_seen)
            max_base_quality_seen = pileup.pile[i].base_q;
        if (pileup.pile[i].map_q < min_map_quality_seen)
            min_map_quality_seen = pileup.pile[i].map_q;
        else if (pileup.pile[i].map_q > max_map_quality_seen)
            max_map_quality_seen = pileup.pile[i].map_q;
    }
    if (debug(2)) {
        if (prev_min_bq > min_base_quality_seen)
            std::cerr << "NL=" << NL << " min_base_quality_seen = " << min_base_quality_seen
                << ":" << uint16_t(min_base_quality_seen) << std::endl;
        if (prev_max_bq > max_base_quality_seen)
            std::cerr << "NL=" << NL << " max_base_quality_seen = " << max_base_quality_seen
                << ":" << uint16_t(max_base_quality_seen) << std::endl;
        if (prev_min_mq > min_map_quality_seen)
            std::cerr << "NL=" << NL << " min_map_quality_seen = " << min_map_quality_seen
                << ":" << uint16_t(min_map_quality_seen) << std::endl;
        if (prev_max_mq > max_map_quality_seen)
            std::cerr << "NL=" << NL << " max_map_quality_seen = " << max_map_quality_seen
                << ":" << uint16_t(max_map_quality_seen) << std::endl;
    }
}


//----------------- printing


void
PileupParser::print(std::ostream& os) const
{
    os << (*this) << std::endl;
}


void
PileupParser::print_lite(std::ostream& os, const std::string sep) const
{
    os << filename << sep << NL << ":" << NF;
    for (size_t i = F_ref; i < F_END; ++i) os << sep << fields[i];
    os << std::endl;
}


std::ostream& 
operator<<(std::ostream&os, const PileupParser& parser)
{
    std::string sep(" ");
    os << parser.filename << ":" << parser.NL << ":" << parser.line << std::endl;
    os << parser.filename << ":" << parser.NL << ":NF=" << parser.NF;
    for (size_t i = 0; i < parser.fields.size(); ++i) os << sep << i << "=" << parser.fields[i];
    return(os);
}


//--------------------------------------------------------
//--------------------------------- class Pileup

// Describes one position in the pileup
//
// ref              : reference sequence name (TODO: make more space-efficient)
// pos              : base position within reference sequence (1-bases)
// cov              : coverage as reported in the pileup (-1 is not set)
// raw_base_call    : pointer to *unparsed* fields[4] member for base calls
// raw_base_quality : pointer to *unparsed* fields[5] member for base quality
// raw_map_quality  : pointer to *unparsed* fields[6] member for mapping quality
// pile             : vector of Stratum, describing each read contribution
// indels           : vector of Indel, describing each declared indel 
// parse_state      : PS_NONE, PS_lite, PS_pile, PS_all (== PS_lite | PS_pile),
//                    tells when the line has been read and the pileup has been
//                    parsed
// min_set_base_quality : set by user
// min_set_map_quality  : set by user, defaults to samtools +33

Pileup::Pileup(uchar_t min_base_qual) 
    : ref(""), pos(0), refbase('\0'), cov(-1), 
      raw_base_call(0), raw_base_quality(0), raw_map_quality(0),
      parse_state(PS_NONE), 
      min_set_base_quality(min_base_qual), min_set_map_quality(33)
{ }


Pileup::~Pileup()
{ }


void
Pileup::reset_pile()
{
    // should never affect anything set by parse_line_lite()
    pile.clear();
    indels.clear();
}


BaseCount
Pileup::base_count()
{
    BaseCount ans;
    for (Pile::const_iterator citer = pile.begin(); citer != pile.end(); ++citer)
        ++ans[citer->base];
    return(ans);
}


bool 
Pileup::set_min_base_quality(uchar_t min_base_q)
{
    if (min_base_q == 0)
        std::cerr << "pileup::set_min_base_quality: argument is 0, so a no-op" << std::endl;
    bool good_quals = true;
    min_set_base_quality = min_base_q;
    for (Pile::iterator iter = pile.begin(); iter != pile.end(); ++iter) {
        if (iter->base_q < min_base_q)
            good_quals = false;
        iter->base_q -= min_base_q;
    }
    return good_quals;
}


bool 
Pileup::set_min_map_quality(uchar_t min_map_q)
{
    if (min_map_q == 0)
        std::cerr << "pileup::set_min_base_quality: argument is 0, so a no-op" << std::endl;
    bool good_quals = true;
    min_set_map_quality = min_map_q;
    for (Pile::iterator iter = pile.begin(); iter != pile.end(); ++iter) {
        if (iter->map_q < min_map_q or iter->read_map_q < min_map_q)
            good_quals = false;
        iter->map_q -= min_map_q;
        if (iter->read_map_q) 
            iter->read_map_q -= min_map_q;
    }
    return good_quals;
}


std::vector<uchar_t>
Pileup::get_map_q(const size_t start, size_t end)
{
    std::vector<uchar_t> ans(pile.size());
    end = std::min(end, pile.size() - 1);
    for (size_t i = start; i <= end; ++i) 
        ans[i] = pile[i].map_q;
    return ans;
}
   

void 
Pileup::print(std::ostream& os) const
{
    os << (*this) << std::endl;
}


void
Pileup::print_pile(std::ostream& os, const size_t start, size_t end, const std::string sep) const
{
    if (pile.size() == 0) {
        os << "NO_PILE"; 
    } else {
        print_pile_stack(os, start, end, true, sep);
    }
}


void
Pileup::print_pile_stack(std::ostream& os, const size_t start, size_t end, 
                         const bool include_indels, const std::string end_stack) const
{
    if (pile.size() == 0) {
        os << "NO_PILE" << end_stack; 
    } else {
        end = std::min(std::max(end, start), pile.size() - 1);
        size_t i;
        for (i = start; i <= end; ++i) {
            os << static_cast<uchar_t>(pile[i].dir == RD_fwd ?  pile[i].base : tolower(pile[i].base));
            if (include_indels && pile[i].indel)
                pile[i].indel->print_compact(os);
        }
        os << end_stack;
        for (i = start; i <= end; ++i) os << pile[i].base_q;
        os << end_stack;
        for (i = start; i <= end; ++i) os << pile[i].map_q;
        os << end_stack;
        os.flush();
    }
}


void 
Pileup::debug_print() 
{ 
    print(std::cerr);
}


void 
Pileup::debug_print_pile(const size_t start, size_t end, bool stack)
{
    if (stack) print_pile_stack(std::cerr, start, end);
    else print_pile(std::cerr, start, end);
}
         


std::ostream& 
operator<<(std::ostream& os, const Pileup& pileup)
{
    std::string sep(" ");
    os << "pileup";
    os << sep << "0x" << pileup.parse_state;
    os << sep << pileup.ref;
    os << sep << pileup.pos;
    os << sep << pileup.refbase;
    os << sep << pileup.cov;
    //os << sep << "[b=" << (uint16_t)min_set_base_quality;
    //os << ",m=" << (uint16_t)min_set_map_quality << "]";
    os << sep;
    pileup.print_pile(os, 0, pileup.pile.size() - 1, sep);
    return(os);
}


//--------------------------------------------------------
//--------------------------------- class Stratum

// One stratum per read per position
//
// base        : the base declared in this stratum
// base_q      : base quality, scale must be set by user (TODO: guess quality, tell user if they're wrong)
// map_q       : mapping quality, samtools convention is Phred+33s (TODO: guess quality, & other mappers?)
// dir         : RD_NONE, RD_fwd, RD_rev (enum typedef in PileupTools namespace)
// read_str    : read structure, RS_NONE, RS_start, RS_end, RS_gap (enum typedef in PileupTools namespace)
// read_map_q  : mapping quality, samtools convention is Phred+33s (TODO: other mappers?)
// indel       : pointer to class Indel instance if there's an indel declared here

Stratum::Stratum() 
    : base('\0'), base_q('\0'), map_q('\0'), dir(RD_NONE), read_str(RS_NONE),
      read_map_q('\0'), indel(0)
{ }

Stratum::~Stratum()
{ }


//--------------------------------------------------------
//--------------------------------- class Indel

// If an indel is declared, an instance of this holds it.  These are
// managed in a vector held in a Pileup, and individual instances of
// Stratum point to the entry in this vector.
//
// type     : IN_NONE, IN_ins, IN_del (enum typedef in PileupTools namespace)
// dir      : RD_NONE, RD_fwd, RD_rev (enum typedef in PileupTools namespace)
// size     : signed, with + = insertion, - = deletion (yes, the sign is redundant to type)
// seq      : the uppercase sequence of the indel
// stratum  : the read stratum in which the Indel was declared
// map_q    : mapping quality of the read declaring the indel
//
// TODO: anything else to note for an indel?

Indel::Indel(const int32_t sz, const std::string& sq, const size_t strat, const uchar_t mq)
    : type(sz > 0 ? IN_ins : IN_del),
      dir(isBaseForward(sq[0]) ? RD_fwd : RD_rev),
      size(sz),
      seq(dir == RD_fwd ? sq : toUpper(sq)),
      stratum(strat),
      map_q(mq)
{ }


Indel::Indel()
    : type(IN_NONE), dir(RD_NONE), size(0), seq(""), stratum(0), map_q(0)
{ }


Indel::~Indel()
{ }


std::ostream& 
operator<<(std::ostream&os, const Indel& indel)
{
    std::string sep(" ");
    os << "(@" << indel.stratum << sep;
    if (indel.type == PileupTools::IN_NONE)
        os << "0!!!)";
    else
        os << (indel.dir == PileupTools::RD_fwd ? "." : ",") << sep << indel.size 
            << sep << indel.seq_qualified() << ")";
    return(os);
}


void
Indel::print(std::ostream& os, const std::string sep) const
{
    os << "indel" << (*this) << std::endl;
}


void
Indel::print_compact(std::ostream& os) const
{
    if (type == IN_NONE)
        os << "(0)";
    else
        os << "(" << seq_qualified() << ")";
}


const std::string
Indel::seq_qualified() const
{
    std::ostringstream qseq;
    qseq << (type > IN_ins ? "+" : "-") << (dir == RD_rev ? toLower(seq) : seq);
    return(qseq.str());
}


} // namespace PileupTools
