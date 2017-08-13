#include <cstdio>
#include <fstream>
#include "io/carp.h"
#include "util/CarpStreamBuf.h"
#include "util/AminoAcidUtil.h"
#include "util/Params.h"
#include "util/FileUtils.h"
#include "util/StringUtils.h"
#include "GeneratePeptides.h"
#include "TideIndexApplication.h"
#include "TideMatchSet.h"
#include "app/tide/modifications.h"
#include "app/tide/records_to_vector-inl.h"
#include <omp.h>
#ifdef _MSC_VER
#include <io.h>
#endif

extern void AddTheoreticalPeaks(const vector<const pb::Protein*>& proteins,
                const string& input_filename, const string& output_filename);
extern void AddMods(HeadedRecordReader* reader, string out_file, string tmpDir,
                const pb::Header& header, const vector<const pb::Protein*>& proteins,
                VariableModTable& var_mod_table);
DECLARE_int32(max_mods);
DECLARE_int32(min_mods);

class WMUPeptidesLocation {
public:
	short thread_number_id;
	int peptide_location;
};

class WMUProtien {
public:
	int protienNumber;
	string protienName;
	string protienSequance;
	vector<WMUPeptidesLocation> peptidesLocationsVector;
};

class WMUPeptide {
public:
	string peptideSequence;
	string decoySequence;
	int protienNumber;
	int peptideNumber;
	float mass;
	int location;
	int decoyProtienNumber;
};

int THREAD_NUM = 1;
vector<WMUProtien> protiensVector;
vector<WMUPeptide> *peptidesVector;

map<string, bool> *generatedPeptidesSet;
set<string> *generatedDecoysSet;

TideIndexApplication::TideIndexApplication() {
}

TideIndexApplication::~TideIndexApplication() {
}

int TideIndexApplication::main(int argc, char** argv) {

        return main(Params::GetString("protein fasta file"),
                        Params::GetString("index name"),
                        StringUtils::Join(vector < string > (argv, argv + argc), ' '));
}

int TideIndexApplication::main(const string& fasta, const string& index,
                string cmd_line) {

		carp(CARP_INFO, "===============================================================");
		carp(CARP_INFO, "=                     Version 2.0                             =");
		carp(CARP_INFO, "=                     May 22, 2017                            =");
		carp(CARP_INFO, "===============================================================");


        carp(CARP_INFO, "Running tide-index...");

        if (cmd_line.empty()) {
                cmd_line = "crux tide-index " + fasta + " " + index;
        }

        // Reroute stderr
        CarpStreamBuf buffer;
        streambuf* old = cerr.rdbuf();
        cerr.rdbuf(&buffer);

        // Get options
        double min_mass = Params::GetDouble("min-mass");
        double max_mass = Params::GetDouble("max-mass");
        int min_length = Params::GetInt("min-length");
        int max_length = Params::GetInt("max-length");
        bool monoisotopic_precursor = Params::GetString("isotopic-mass")
                        != "average";
        FLAGS_max_mods = Params::GetInt("max-mods");
        FLAGS_min_mods = Params::GetInt("min-mods");
        bool allowDups = Params::GetBool("allow-dups");
        if (FLAGS_min_mods > FLAGS_max_mods) {
                carp(CARP_FATAL,
                                "The value for 'min-mods' cannot be greater than the value "
                                                "for 'max-mods'");
        }
        MASS_TYPE_T mass_type = (monoisotopic_precursor) ? MONO : AVERAGE;
        int missed_cleavages = Params::GetInt("missed-cleavages");
        DIGEST_T digestion = get_digest_type_parameter("digestion");
        ENZYME_T enzyme_t = get_enzyme_type_parameter("enzyme");
        const char* enzymePtr = enzyme_type_to_string(enzyme_t);
        string enzyme(enzymePtr);
        if ((enzyme != "no-enzyme")
                        && (digestion != FULL_DIGEST && digestion != PARTIAL_DIGEST)) {
                carp(CARP_FATAL,
                                "'digestion' must be 'full-digest' or 'partial-digest'");
        }

        VariableModTable var_mod_table;
        var_mod_table.ClearTables();
        //parse regular amino acid modifications
        string mods_spec = Params::GetString("mods-spec");
        carp(CARP_DEBUG, "mods_spec='%s'", mods_spec.c_str());
        if (!var_mod_table.Parse(mods_spec.c_str())) {
                carp(CARP_FATAL, "Error parsing mods");
        }
        //parse terminal modifications
        mods_spec = Params::GetString("cterm-peptide-mods-spec");
        if (!mods_spec.empty() && !var_mod_table.Parse(mods_spec.c_str(), CTPEP)) {
                carp(CARP_FATAL, "Error parsing c-terminal peptide mods");
        }
        mods_spec = Params::GetString("nterm-peptide-mods-spec");
        if (!mods_spec.empty() && !var_mod_table.Parse(mods_spec.c_str(), NTPEP)) {
                carp(CARP_FATAL, "Error parsing n-terminal peptide mods");
        }
        mods_spec = Params::GetString("cterm-protein-mods-spec");
        if (!mods_spec.empty() && !var_mod_table.Parse(mods_spec.c_str(), CTPRO)) {
                carp(CARP_FATAL, "Error parsing c-terminal protein mods");
        }
        mods_spec = Params::GetString("nterm-protein-mods-spec");
        if (!mods_spec.empty() && !var_mod_table.Parse(mods_spec.c_str(), NTPRO)) {
                carp(CARP_FATAL, "Error parsing n-terminal protein mods");
        }

        var_mod_table.SerializeUniqueDeltas();

        if (!MassConstants::Init(var_mod_table.ParsedModTable(),
                        var_mod_table.ParsedNtpepModTable(),
                        var_mod_table.ParsedCtpepModTable(), 0, 0)) {
                carp(CARP_FATAL, "Error in MassConstants::Init");
        }

        DECOY_TYPE_T decoy_type = get_tide_decoy_type_parameter("decoy-format");
        string decoyPrefix = Params::GetString("decoy-prefix");

        // Set up output paths
        bool overwrite = Params::GetBool("overwrite");

        if (!FileUtils::Exists(fasta)) {
                carp(CARP_FATAL, "Fasta file %s does not exist", fasta.c_str());
        }

        string out_proteins = FileUtils::Join(index, "protix");
        string out_peptides = FileUtils::Join(index, "pepix");
        string out_aux = FileUtils::Join(index, "auxlocs");
        string modless_peptides = out_peptides + ".nomods.tmp";
        string peakless_peptides = out_peptides + ".nopeaks.tmp";
        ofstream* out_target_list = NULL;
        ofstream* out_decoy_list = NULL;
        if (Params::GetBool("peptide-list")) {
                out_target_list = create_stream_in_path(
                                make_file_path("tide-index.peptides.target.txt").c_str(), NULL,
                                overwrite);
                if (decoy_type != NO_DECOYS) {
                        out_decoy_list = create_stream_in_path(
                                        make_file_path("tide-index.peptides.decoy.txt").c_str(),
                                        NULL, overwrite);
                }
        }
        ofstream* out_decoy_fasta =
                        GeneratePeptides::canGenerateDecoyProteins() ?
                                        create_stream_in_path(
                                                        make_file_path("tide-index.decoy.fasta").c_str(),
                                                        NULL, overwrite) :
                                        NULL;

        if (create_output_directory(index.c_str(), overwrite) != 0) {
                carp(CARP_FATAL, "Error creating index directory");
        } else if (FileUtils::Exists(out_proteins)
                        || FileUtils::Exists(out_peptides) || FileUtils::Exists(out_aux)) {
                if (overwrite) {
                        carp(CARP_DEBUG, "Cleaning old index file(s)");
                        FileUtils::Remove(out_proteins);
                        FileUtils::Remove(out_peptides);
                        FileUtils::Remove(out_aux);
                        FileUtils::Remove(modless_peptides);
                        FileUtils::Remove(peakless_peptides);
                } else {
                        carp(CARP_FATAL,
                                        "Index file(s) already exist, use --overwrite T or a "
                                                        "different index name");
                }
        }

        // Start tide-index
        carp(CARP_INFO, "Reading %s and computing unmodified peptides...",
                        fasta.c_str());
        pb::Header proteinPbHeader;

        vector<string*> proteinSequences;

         struct timeval after, before;
          gettimeofday(&before, 0);
          THREAD_NUM = Params::GetInt("num-threads");
          if (THREAD_NUM == 0) {
        	  THREAD_NUM = 1;
          }
          peptidesVector = new vector<WMUPeptide>[THREAD_NUM];



          generatedPeptidesSet = new map<string, bool>[THREAD_NUM];
          generatedDecoysSet = new set<string>[THREAD_NUM];


//          fastaToPb_OLD(cmd_line, enzyme_t, digestion, missed_cleavages, min_mass,
//                        max_mass, min_length, max_length, allowDups, mass_type, decoy_type,
//                        fasta, out_proteins, proteinPbHeader, outPeptideHeap, proteinSequences,
//                        out_decoy_fasta);
          fastaToPb(cmd_line, enzyme_t, digestion, missed_cleavages, min_mass,
                                  max_mass, min_length, max_length, allowDups, mass_type, decoy_type,
                                  fasta, out_proteins, proteinPbHeader, proteinSequences,
                                  out_decoy_fasta);



          int counter = 0;

          gettimeofday(&after, 0);
          double elapsed = (after.tv_sec - before.tv_sec) + (after.tv_usec - before.tv_usec) / 1000000.0;

          carp(CARP_INFO, "Elapsed Time inside fastaToPb %f", elapsed);

        pb::Header header_with_mods;

        // Set up peptides header
        pb::Header_PeptidesHeader& pep_header =
                        *(header_with_mods.mutable_peptides_header());
        pep_header.Clear();
        pep_header.set_min_mass(min_mass);
        pep_header.set_max_mass(max_mass);
        pep_header.set_min_length(min_length);
        pep_header.set_max_length(max_length);
        pep_header.set_monoisotopic_precursor(monoisotopic_precursor);
        pep_header.set_enzyme(enzyme);
        if (enzyme != "no-enzyme") {
                pep_header.set_full_digestion(digestion == FULL_DIGEST);
                pep_header.set_max_missed_cleavages(missed_cleavages);
        }
        pep_header.mutable_mods()->CopyFrom(*(var_mod_table.ParsedModTable()));
        pep_header.mutable_nterm_mods()->CopyFrom(
                        *(var_mod_table.ParsedNtpepModTable()));
        pep_header.mutable_cterm_mods()->CopyFrom(
                        *(var_mod_table.ParsedCtpepModTable()));

        header_with_mods.set_file_type(pb::Header::PEPTIDES);
        header_with_mods.set_command_line(cmd_line);
        pb::Header_Source* source = header_with_mods.add_source();
        source->mutable_header()->CopyFrom(proteinPbHeader);
        source->set_filename(AbsPath(out_proteins));

        pb::Header header_no_mods;
        header_no_mods.CopyFrom(header_with_mods);
        pb::ModTable* del =
                        header_no_mods.mutable_peptides_header()->mutable_mods();
        del->mutable_variable_mod()->Clear();
        del->mutable_unique_deltas()->Clear();

        bool need_mods = var_mod_table.Unique_delta_size() > 0;

        string basic_peptides = need_mods ? modless_peptides : peakless_peptides;
        carp(CARP_DETAILED_DEBUG, "basic_peptides=%s", basic_peptides.c_str());

        writePeptidesAndAuxLocs(basic_peptides, out_aux,
                        header_no_mods);
        // Do some clean up
        for (vector<string*>::iterator i = proteinSequences.begin();
                        i != proteinSequences.end(); ++i) {
                delete *i;
        }
        vector<TideIndexPeptide>().swap(outPeptideHeap);
        ProteinVec proteins;
        if (!ReadRecordsToVector<pb::Protein>(&proteins, out_proteins)) {
                carp(CARP_FATAL, "Error reading proteins file");
        }

        if (need_mods) {
                carp(CARP_INFO, "Computing modified peptides...");
                HeadedRecordReader reader(modless_peptides, NULL, 1024 << 10); // 1024kb buffer
                AddMods(&reader, peakless_peptides, Params::GetString("temp-dir"),
                                header_with_mods, proteins, var_mod_table);
        }

        if (out_target_list) {
                // Write peptide lists
                carp(CARP_INFO, "Writing peptide lists...");

                // This set holds target peptide strings
                set<string> targetPepStrs;
                // This vector holds decoy peptide strings and masses
                vector < pair<string, double> > decoyPepStrs;
                // Read peptides protocol buffer file
                vector<const pb::AuxLocation*> locations;
                if (!ReadRecordsToVector<pb::AuxLocation>(&locations, out_aux)) {
                        carp(CARP_FATAL, "Error reading auxlocs file");
                }
                // Iterate over all protocol buffer peptides
                unsigned int writeCountTargets = 0, writeCountDecoys = 0;
                HeadedRecordReader reader(peakless_peptides, NULL);
                while (!reader.Done()) {
                        pb::Peptide* protobuf = new pb::Peptide;
                        reader.Read(protobuf);
                        pb::Peptide* peptide = protobuf;
                        bool writeTarget = true;
                        bool writeDecoy = false;

                        if (out_decoy_list) {
                                writeDecoy = peptide->is_decoy();
                                writeTarget = !writeDecoy;
                        }
                        string pep_str = getModifiedPeptideSeq(peptide, &proteins);

                        if (writeTarget) {
                                // This is a target, output it
                                targetPepStrs.insert(pep_str);
                                *out_target_list << pep_str << '\t' << peptide->mass() << endl;
                                ++writeCountTargets;
                        }
                        if (writeDecoy) {
                                // This is a decoy, save it to output later
                                decoyPepStrs.push_back(make_pair(pep_str, peptide->mass()));
                                ++writeCountDecoys;
                        }
                        delete peptide;
                }

                // Iterate over saved decoys and output them
                for (vector<pair<string, double> >::iterator i = decoyPepStrs.begin();
                                i != decoyPepStrs.end(); ++i) {
                        *out_decoy_list << i->first << '\t' << i->second;
                        if (targetPepStrs.find(i->first) != targetPepStrs.end()) {
                                *out_decoy_list << "\t*";
                        }
                        *out_decoy_list << endl;
                }

                // Close and clean up streams
                if (out_decoy_list) {
                        out_decoy_list->close();
                        delete out_decoy_list;
                }
                out_target_list->close();
                delete out_target_list;

                carp(CARP_DEBUG, "Wrote %d targets and %d decoys to peptide list",
                                writeCountTargets, writeCountDecoys);
        }

        carp(CARP_INFO, "Precomputing theoretical spectra...");
        AddTheoreticalPeaks(proteins, peakless_peptides, out_peptides);

        // Clean up
        for (vector<const pb::Protein*>::iterator i = proteins.begin();
                        i != proteins.end(); ++i) {
                delete *i;
        }

        // clean up out_decoy_fasta
        if (out_decoy_fasta) {
                delete out_decoy_fasta;
        }

        // Recover stderr
        cerr.rdbuf(old);
        FileUtils::Remove(modless_peptides);
        FileUtils::Remove(peakless_peptides);

        return 0;
}

string TideIndexApplication::getName() const {
        return "tide-index";
}

string TideIndexApplication::getDescription() const {
        return "[[nohtml:Create an index for all peptides in a fasta file, for use in "
                        "subsequent calls to tide-search.]]"
                        "[[html:<p>Tide is a tool for identifying peptides from tandem mass "
                        "spectra. It is an independent reimplementation of the SEQUEST<sup>&reg;"
                        "</sup> algorithm, which assigns peptides to spectra by comparing the "
                        "observed spectra to a catalog of theoretical spectra derived from a "
                        "database of known proteins. Tide's primary advantage is its speed. Our "
                        "published paper provides more detail on how Tide works. If you use Tide "
                        "in your research, please cite:</p><blockquote>Benjamin J. Diament and "
                        "William Stafford Noble. &quot;<a href=\""
                        "http://dx.doi.org/10.1021/pr101196n\">Faster SEQUEST Searching for "
                        "Peptide Identification from Tandem Mass Spectra.</a>&quot; <em>Journal of "
                        "Proteome Research</em>. 10(9):3871-9, 2011.</blockquote><p>The <code>"
                        "tide-index</code> command performs an optional pre-processing step on the "
                        "protein database, converting it to a binary format suitable for input to "
                        "the <code>tide-search</code> command.</p><p>Tide considers only the "
                        "standard set of 20 amino acids. Peptides containing non-amino acid "
                        "alphanumeric characters (BJOUXZ) are skipped. Non-alphanumeric characters "
                        "are ignored completely.</p>]]";
}

vector<string> TideIndexApplication::getArgs() const {
        string arr[] = { "protein fasta file", "index name" };
        return vector < string > (arr, arr + sizeof(arr) / sizeof(string));
}

vector<string> TideIndexApplication::getOptions() const {
        string arr[] = { "decoy-format", "keep-terminal-aminos", "decoy-prefix",
                        "enzyme", "custom-enzyme", "digestion", "missed-cleavages",
                        "max-length", "max-mass", "min-length", "min-mass", "isotopic-mass",
                        "mods-spec", "cterm-peptide-mods-spec", "nterm-peptide-mods-spec",
                        "max-mods", "min-mods", "output-dir", "overwrite", "peptide-list",
                        "parameter-file", "seed", "clip-nterm-methionine", "verbosity",
                        "allow-dups", "temp-dir" };
        return vector < string > (arr, arr + sizeof(arr) / sizeof(string));
}

vector<pair<string, string> > TideIndexApplication::getOutputs() const {
        vector < pair<string, string> > outputs;
        outputs.push_back(
                        make_pair("index",
                                        "A binary index, using the name specified on the command line."));
        outputs.push_back(
                        make_pair("tide-index.params.txt",
                                        "a file containing the name and value of all parameters/options for the "
                                                        "current operation. Not all parameters in the file may have been used in "
                                                        "the operation. The resulting file can be used with the --parameter-file "
                                                        "option for other crux programs."));
        outputs.push_back(
                        make_pair("tide-index.log.txt",
                                        "a log file containing a copy of all messages that were printed to the "
                                                        "screen during execution."));
        return outputs;
}

bool TideIndexApplication::needsOutputDirectory() const {
        return true;
}

COMMAND_T TideIndexApplication::getCommand() const {
        return TIDE_INDEX_COMMAND;
}

//#define THREAD_NUM = PARAM
#define REAL_PEPTIDE true
#define DECOY_PEPTIDE false

/*
static bool makeDecoy(
  const std::string& seq, ///< sequence to make decoy from
  const std::set<std::string>& targetSeqs,  ///< targets to check against
  const std::set<std::string>& decoySeqs,  ///< decoys to check against
  bool shuffle, ///< shuffle (if false, reverse)
  std::string& decoyOut ///< string to store decoy
);
*/
bool shufflePeptide(
  string& seq, ///< Peptide sequence to shuffle
  unsigned int maxShuffleAttempts ///< Maximum number of shuffle attempts
) {
  switch (seq.length()) {
  case 0:
  case 1:
    return false;
  case 2: {
    char second = seq[1];
    seq[1] = seq[0];
    seq[0] = second;
    return true;
  }
  default:
    string originalSeq(seq);
    for (int i = 0; i < maxShuffleAttempts; i++) {
      random_shuffle(seq.begin(), seq.end(), myrandom_limit);
      if (seq != originalSeq) {
        return true;
      }
    }
    return false;
  }
}

/**
 * Reverses the peptide sequence.
 * Returns false if no different sequence was generated
 */
bool reversePeptide(
  string& seq ///< Peptide sequence to reverse
) {
  string originalSeq(seq);
  reverse(seq.begin(), seq.end());
  return seq != originalSeq;
}
/*
bool wmuMakeDecoy(const string& seq, ///< sequence to make decoy from
		const map<string, bool>& targetSeqs, ///< targets to check against
		const set<string>& decoySeqs, ///< decoys to check against
		bool shuffle, ///< shuffle (if false, reverse)
		string& decoyOut ///< string to store decoy
		) {

	  const string keepTerminal = Params::GetString("keep-terminal-aminos");
	  string decoyPre, decoyPost;
	  if (keepTerminal == "N") {
	    if (seq.length() <= 2) {
	      decoyOut = seq;
	      return false;
	    }
	    decoyPre = seq[0];
	    decoyOut = seq.substr(1);
	  } else if (keepTerminal == "C") {
	    if (seq.length() <= 2) {
	      decoyOut = seq;
	      return false;
	    }
	    decoyPost = seq[seq.length() - 1];
	    decoyOut = seq.substr(0, seq.length() - 1);
	  } else if (keepTerminal == "NC") {
	    if (seq.length() <= 3) {
	      decoyOut = seq;
	      return false;
	    }
	    decoyPre = seq[0];
	    decoyPost = seq[seq.length() - 1];
	    decoyOut = seq.substr(1, seq.length() - 2);
	  } else {
	    decoyOut = seq;
	    if (seq.length() <= 1) {
	      return false;
	    }
	  }

	  if (!shuffle) {
	    // Reverse
	    if (reversePeptide(decoyOut)) {
	      // Re-add n/c
	      string decoyCheck = decoyPre + decoyOut + decoyPost;
	      // Check in sets
	      if (targetSeqs.find(decoyCheck) == targetSeqs.end() &&
	          decoySeqs.find(decoyCheck) == decoySeqs.end()) {
	        decoyOut = decoyCheck;
	        return true;
	      }
	    }
	    carp(CARP_DEBUG, "Failed reversing %s, shuffling", seq.c_str());
	  }

	  // Shuffle
	  if (shufflePeptide(decoyOut, 6)) {
	    // Re-add n/c
	    string decoyCheck = decoyPre + decoyOut + decoyPost;
	    // Check in sets
	    if (targetSeqs.find(decoyCheck) == targetSeqs.end() &&
	        decoySeqs.find(decoyCheck) == decoySeqs.end()) {
	      decoyOut = decoyCheck;
	      return true;
	    }
	  }

	  decoyOut = seq;
	  return false;
}*/
bool wmuMakeDecoy(const string& seq, ///< sequence to make decoy from
                const map<string, bool>& targetSeqs, ///< targets to check against
                const set<string>& decoySeqs, ///< decoys to check against
                bool shuffle, ///< shuffle (if false, reverse)
                string& decoyOut ///< string to store decoy
               ) {

		const string keepTerminal = Params::GetString("keep-terminal-aminos");
        if (seq.empty() || seq.length() <= 2) {
                return false;
        }

        string decoyPre, decoyPost;
        if (keepTerminal == "N") {
                decoyPre = seq[0];
                decoyOut = seq.substr(1);
        } else if (keepTerminal == "C") {
                decoyPost = seq[seq.length() - 1];
                decoyOut = seq.substr(0, seq.length() - 1);
        } else if (keepTerminal == "NC") {
                decoyPre = seq[0];
                decoyPost = seq[seq.length() - 1];
                decoyOut = seq.substr(1, seq.length() - 2);
                //decoyOut="LPETIK";
        } else {
                decoyOut = seq;
        }

        if (shuffle) {
                bool isShuffled = false;
                //string originalSeq(decoyOut);
				for (int i = 0; i < decoyOut.length(); i += 2) {
						swap(decoyOut[i], decoyOut[i + 1]);
				}
            
                int random1 = rand()%(decoyOut.length() - 0 + 1) + 0;
                int random2 = rand()%(decoyOut.length() - 0 + 1) + 0;
                
                swap(decoyOut[random1], decoyOut[random2]);
				decoyOut = decoyPre + decoyOut + decoyPost;
                if (targetSeqs.find(decoyOut) == targetSeqs.end()
                                && decoySeqs.find(decoyOut) == decoySeqs.end()) {
                        return true;
                } else {
                        return false;
                }
        } else {
                //reverse(seq.begin(), seq.end());
        }
        return true;
}

void TideIndexApplication::fastaToPb(
                const string& commandLine,
                const ENZYME_T enzyme,
                const DIGEST_T digestion,
                int missedCleavages,
                FLOAT_T minMass,
                FLOAT_T maxMass,
                int minLength,
                int maxLength,
                bool allowDups,
                MASS_TYPE_T massType,
                DECOY_TYPE_T decoyType,
                const string& fasta,
                const string& proteinPbFile,
                pb::Header& outProteinPbHeader,
                //vector<TideIndexPeptide>& outPeptideHeap,
                vector<string*>& outProteinSequences,
                ofstream* decoyFasta
                ) {

	string decoyPrefix = Params::GetString("decoy-prefix");
	string keepTerminal = Params::GetString("keep-terminal-aminos");
	outProteinPbHeader.Clear();
	outProteinPbHeader.set_file_type(pb::Header::RAW_PROTEINS);
	outProteinPbHeader.set_command_line(commandLine);
	pb::Header_Source* headerSource = outProteinPbHeader.add_source();
	HeadedRecordWriter proteinWriter(proteinPbFile, outProteinPbHeader);
	pb::Protein pbProtein;
	outPeptideHeap.clear();

	unsigned int targetsGenerated[THREAD_NUM];
	unsigned int decoysGenerated[THREAD_NUM];
	unsigned int invalidPepCnt[THREAD_NUM];
	unsigned int failedDecoyCnt[THREAD_NUM];

	unsigned int totalTargetsGenerated = 0;
	unsigned int totalDecoysGenerated = 0;
	unsigned int totalInvalidPepCnt = 0;
	unsigned int totalFailedDecoyCnt = 0;

	typedef GeneratePeptides::CleavedPeptide PeptideInfo;
	ifstream fastaStream(fasta.c_str(), ifstream::in);
	string line;

	double proteinNumberSequence = 0;
	for (int i = 0; i < THREAD_NUM; i++) {
		targetsGenerated[i] = 0;
		decoysGenerated[i] = 0;
		invalidPepCnt[i] = 0;
		failedDecoyCnt[i] = 0;
	}
    
	while (true) {
		WMUProtien protien;
		protien.protienNumber = protiensVector.size();
		if (!GeneratePeptides::getNextProtein(fastaStream, &(protien.protienName), &(protien.protienSequance))) {
			break;
		}
		protiensVector.push_back(protien);
	}


	/*
	while (fastaStream.good() && !fastaStream.eof()) {
		getline(fastaStream, line);
		line = StringUtils::Trim(line);
		if (StringUtils::StartsWith(line, ">")) {
			WMUProtien protien;

			protien.protienNumber = protiensVector.size();
			bool idStart = false;
			string::const_iterator begin = line.begin() + 1;
			string::const_iterator end = line.end();

			for (string::const_iterator i = begin; i != line.end(); i++) {
				bool space = isspace(*i);
				if (!idStart && !space) {
					begin = i;
					idStart = true;
				} else if (idStart && space) {
					end = i;
					break;
				}
			}

			protien.protienName = string(begin, end);
			while (fastaStream.good()) {
				if (fastaStream.eof() || fastaStream.peek() == '>') {
						break;
				}
				getline(fastaStream, line);
				protien.protienSequance.append(string(line));
			}

			protiensVector.push_back(protien);
		}
	}*/

#pragma omp parallel num_threads(THREAD_NUM)
	{

		int thread_number = omp_get_thread_num();
		int peptideNumber = 0;
		vector<PeptideInfo> cleavedPeptides;

#pragma omp for
		for (vector<WMUProtien>::iterator it = protiensVector.begin(); it < protiensVector.end(); it++) {

			cleavedPeptides = GeneratePeptides::cleaveProtein(it->protienSequance, enzyme, digestion, missedCleavages, minLength, maxLength);

			//cleavedPeptides = GeneratePeptides::cleaveProtein(*proteinSequence   , enzyme, digestion, missedCleavages, minLength, maxLength);
			for (vector<PeptideInfo>::iterator i = cleavedPeptides.begin(); i < cleavedPeptides.end();) {

				FLOAT_T pepMass = calcPepMassTide(i->Sequence(), massType);
				if (pepMass < 0.0) {
					// Sequence contained some invalid character
					carp(CARP_INFO, "Ignoring invalid sequence <%s>", i->Sequence().c_str());
					++invalidPepCnt[thread_number];
					i = cleavedPeptides.erase(i);
					continue;
				} else if (pepMass < minMass || pepMass > maxMass) {
					// Skip to next peptide if not in mass range
					++i;
					continue;
				}

				// Add target to heap
				WMUPeptide peptide;

				peptide.location = i->Position();
				peptide.mass = pepMass;
				peptide.peptideNumber = peptideNumber++;
				peptide.peptideSequence = string(i->Sequence());
				peptide.decoyProtienNumber = 0;
				peptide.protienNumber = it->protienNumber;
				WMUPeptidesLocation peptideLocation;
				int thread_number_id = (int) pepMass % THREAD_NUM;
				peptideLocation.thread_number_id = thread_number_id;

#pragma omp critical
				{
					peptideLocation.peptide_location = peptidesVector[thread_number_id].size();
					peptidesVector[thread_number_id].push_back(peptide);
				}

				it->peptidesLocationsVector.push_back(peptideLocation);
				++targetsGenerated[thread_number];
				++i;
			}
			cleavedPeptides.clear();
		}
	}
    
#pragma omp parallel num_threads(THREAD_NUM)
	{
		int thread_number = omp_get_thread_num();
		string decoySequence;

		for (vector<WMUPeptide>::iterator it = peptidesVector[thread_number].begin(); it < peptidesVector[thread_number].end(); it++) {
			generatedPeptidesSet[thread_number].insert(std::pair<string, string *>(it->peptideSequence, false));
		}

		// generate decoys
		if (decoyType == PROTEIN_REVERSE_DECOYS) {
			// Please use the sequential code 
		} else if (!allowDups) {

			for (vector<WMUPeptide>::iterator it = peptidesVector[thread_number].begin(); it < peptidesVector[thread_number].end(); it++) {

				decoySequence = it->peptideSequence;

				map<string, bool>::iterator peptideSequanceIterator = generatedPeptidesSet[thread_number].find(it->peptideSequence);

				if (peptideSequanceIterator != generatedPeptidesSet[thread_number].end() && peptideSequanceIterator->second == false) {

					bool makingDecoyResult = wmuMakeDecoy(it->peptideSequence, generatedPeptidesSet[thread_number], generatedDecoysSet[thread_number], decoyType == PEPTIDE_SHUFFLE_DECOYS, decoySequence);

					if (makingDecoyResult) {

						it->decoySequence = decoySequence;
						generatedDecoysSet[thread_number].insert(decoySequence);
						decoysGenerated[thread_number]++;
						peptideSequanceIterator->second = true;
					} else {
						failedDecoyCnt[thread_number]++;
					}
				}
			}
		} else {
			// Please use the sequential code 
		}
	}

	// Wrap up
//#pragma omp parallel num_threads(2)
	{
	    int thread_number = omp_get_thread_num();
		int counter = 0;
		//if (thread_number == 0)
		{
			for (vector<WMUProtien>::iterator it = protiensVector.begin(); it < protiensVector.end(); it++) {
				getPbProtein(it->protienNumber, it->protienName, it->protienSequance, pbProtein);
				proteinWriter.Write(&pbProtein);
			}

			for (vector<WMUProtien>::iterator it = protiensVector.begin(); it < protiensVector.end(); it++) {
				for (vector<WMUPeptidesLocation>::iterator j = it->peptidesLocationsVector.begin(); j != it->peptidesLocationsVector.end(); j++) {

					if(peptidesVector[j->thread_number_id][j->peptide_location].decoySequence.length() != 0){
						int decoyProtienNumber = protiensVector.size() + counter++;
						getDecoyPbProtein(decoyProtienNumber,
							  ProteinInfo(it->protienName, &(it->protienSequance)),
							  peptidesVector[j->thread_number_id][j->peptide_location].decoySequence,
							  peptidesVector[j->thread_number_id][j->peptide_location].location,
							  pbProtein);

						peptidesVector[j->thread_number_id][j->peptide_location].decoyProtienNumber = decoyProtienNumber;
						proteinWriter.Write(&pbProtein);

						(*decoyFasta) << ">" << decoyPrefix << it->protienName << endl << peptidesVector[j->thread_number_id][j->peptide_location].decoySequence << endl;
					}
				}
			}
		}
		//else if (thread_number == 1)
		{
			
			for (int i = 0; i < THREAD_NUM; i++) {
				for (vector<WMUPeptide>::iterator it = peptidesVector[i].begin(); it < peptidesVector[i].end(); it++) {
					WMUProtien p = protiensVector[it->protienNumber];
					TideIndexPeptide pepTarget(it->mass, it->peptideSequence.length(), &(it->peptideSequence), it->protienNumber, it->location, false);
					outPeptideHeap.push_back(pepTarget);
					push_heap(outPeptideHeap.begin(), outPeptideHeap.end(), greater<TideIndexPeptide>());
				}
			}

			for (int i = 0; i < THREAD_NUM; i++) {

				totalTargetsGenerated += targetsGenerated[i];
				totalDecoysGenerated += decoysGenerated[i];
				totalInvalidPepCnt += invalidPepCnt[i];
				totalFailedDecoyCnt += failedDecoyCnt[i];

				for (vector<WMUPeptide>::iterator it = peptidesVector[i].begin(); it < peptidesVector[i].end(); it++) {
					if (it->decoySequence.length() != 0) {
						WMUProtien p = protiensVector[it->protienNumber];
						
						TideIndexPeptide pepDecoy(it->mass, it->decoySequence.length(), &(it->decoySequence), it->decoyProtienNumber, (it->location > 0) ? 1 : 0, true);
						outPeptideHeap.push_back(pepDecoy);
						push_heap(outPeptideHeap.begin(), outPeptideHeap.end(), greater<TideIndexPeptide>());
					}
				}
			}
		}
	}

	if (totalInvalidPepCnt > 0) {
		carp(CARP_INFO, "Ignoring %d peptide sequences containing unrecognized characters", totalInvalidPepCnt);
	}
	if (totalFailedDecoyCnt > 0) {
		carp(CARP_INFO, "Failed to generate decoys for %d low complexity peptides", totalFailedDecoyCnt);
	}
	carp(CARP_INFO, "FASTA produced %d targets and %d decoys", totalTargetsGenerated, totalDecoysGenerated);
}

void TideIndexApplication::fastaToPb_OLD(const string& commandLine,
                const ENZYME_T enzyme, const DIGEST_T digestion, int missedCleavages,
                FLOAT_T minMass, FLOAT_T maxMass, int minLength, int maxLength,
                bool allowDups, MASS_TYPE_T massType, DECOY_TYPE_T decoyType,
                const string& fasta, const string& proteinPbFile,
                pb::Header& outProteinPbHeader,
                vector<TideIndexPeptide>& outPeptideHeap,
                vector<string*>& outProteinSequences, ofstream* decoyFasta) {
        typedef GeneratePeptides::CleavedPeptide PeptideInfo;
        string decoyPrefix = Params::GetString("decoy-prefix");
        outProteinPbHeader.Clear();
        outProteinPbHeader.set_file_type(pb::Header::RAW_PROTEINS);
        outProteinPbHeader.set_command_line(commandLine);
        pb::Header_Source* headerSource = outProteinPbHeader.add_source();
        headerSource->set_filename(AbsPath(fasta));
        headerSource->set_filetype("fasta");
        unsigned int invalidPepCnt = 0;
        unsigned int failedDecoyCnt = 0;

        outPeptideHeap.clear();
        outProteinSequences.clear();

        HeadedRecordWriter proteinWriter(proteinPbFile, outProteinPbHeader);
        ifstream fastaStream(fasta.c_str(), ifstream::in);
        pb::Protein pbProtein;
        string proteinName;
        string* proteinSequence = new string;
        int curProtein = -1;
        vector < pair<ProteinInfo, vector<PeptideInfo> > > cleavedPeptideInfo;
        set<string> setTargets, setDecoys;
        map<const string*, TargetInfo> targetInfo;

        // Iterate over all proteins in FASTA file
        unsigned int targetsGenerated = 0, decoysGenerated = 0;
        while (GeneratePeptides::getNextProtein(fastaStream, &proteinName,
                        proteinSequence)) {
        	
                outProteinSequences.push_back(proteinSequence);
                cleavedPeptideInfo.push_back(
                                make_pair(ProteinInfo(proteinName, proteinSequence),
                                                vector<PeptideInfo>()));
                const ProteinInfo& proteinInfo = cleavedPeptideInfo.back().first;
                vector<PeptideInfo>& cleavedPeptides = cleavedPeptideInfo.back().second;
                // Write pb::Protein
                getPbProtein(++curProtein, proteinName, *proteinSequence, pbProtein);
                proteinWriter.Write(&pbProtein);
                cleavedPeptides = GeneratePeptides::cleaveProtein(*proteinSequence,
                                enzyme, digestion, missedCleavages, minLength, maxLength);
                // Iterate over all generated peptides for this protein
                for (vector<PeptideInfo>::iterator i = cleavedPeptides.begin();
                                i != cleavedPeptides.end();) {


                        FLOAT_T pepMass = calcPepMassTide(i->Sequence(), massType);
                        
                        if (pepMass < 0.0) {
                                // Sequence contained some invalid character
                                carp(CARP_DEBUG, "Ignoring invalid sequence <%s>",
                                                i->Sequence().c_str());
                                ++invalidPepCnt;
                                i = cleavedPeptides.erase(i);
                                continue;
                        } else if (pepMass < minMass || pepMass > maxMass) {
                                // Skip to next peptide if not in mass range
                                ++i;
                                continue;
                        }
                        // Add target to heap
                        TideIndexPeptide pepTarget(pepMass, i->Length(), proteinSequence,
                                        curProtein, i->Position(), false);
                        outPeptideHeap.push_back(pepTarget);
                        push_heap(outPeptideHeap.begin(), outPeptideHeap.end(),
                                        greater<TideIndexPeptide>());
                        if (!allowDups && decoyType != NO_DECOYS) {
                                const string* setTarget =
                                                &*(setTargets.insert(i->Sequence()).first);
                                targetInfo.insert(
                                                make_pair(setTarget,
                                                                TargetInfo(proteinInfo, i->Position(),
                                                                                pepMass)));
                        }
                        ++targetsGenerated;
                        ++i;
                }
                proteinSequence = new string;
        }
        delete proteinSequence;
        if (targetsGenerated == 0) {
                carp(CARP_FATAL,
                                "No target sequences generated.  Is \'%s\' a FASTA file?",
                                fasta.c_str());
        }

        // Generate decoys
        map<const string, const string*> targetToDecoy;
        if (decoyType == PROTEIN_REVERSE_DECOYS) {
                if (decoyFasta) {
                        carp(CARP_INFO, "Writing reverse-protein fasta and decoys...");
                }
                for (vector<pair<ProteinInfo, vector<PeptideInfo> > >::const_iterator i =
                                cleavedPeptideInfo.begin(); i != cleavedPeptideInfo.end();
                                ++i) {
                        string decoyProtein = *(i->first.sequence);
                        reverse(decoyProtein.begin(), decoyProtein.end());
                        if (decoyFasta) {
                                (*decoyFasta) << ">" << decoyPrefix << i->first.name << endl
                                                << decoyProtein << endl;
                        }
                        vector<PeptideInfo> cleavedReverse =
                                        GeneratePeptides::cleaveProtein(decoyProtein, enzyme,
                                                        digestion, missedCleavages, minLength, maxLength);
                        // Iterate over all generated peptides for this protein
                        for (vector<PeptideInfo>::iterator j = cleavedReverse.begin();
                                        j != cleavedReverse.end(); ++j) {
                                FLOAT_T pepMass = calcPepMassTide(j->Sequence(), massType);
                                if (pepMass < 0.0) {
                                        // Sequence contained some invalid character
                                        carp(CARP_DEBUG,
                                                        "Ignoring invalid sequence in decoy fasta <%s>",
                                                        j->Sequence().c_str());
                                        ++invalidPepCnt;
                                        continue;
                                } else if (pepMass < minMass || pepMass > maxMass) {
                                        // Skip to next peptide if not in mass range
                                        continue;
                                } else if (!allowDups
                                                && setTargets.find(j->Sequence()) != setTargets.end()) {
                                        // Sequence already exists as a target
                                        continue;
                                }
                                string* decoySequence = new string(j->Sequence());
                                outProteinSequences.push_back(decoySequence);

                                // Write pb::Protein
                                getDecoyPbProtein(++curProtein,
                                                ProteinInfo(i->first.name, &decoyProtein),
                                                *decoySequence, j->Position(), pbProtein);
                                proteinWriter.Write(&pbProtein);
                                // Add decoy to heap
                                TideIndexPeptide pepDecoy(pepMass, j->Length(), decoySequence,
                                                curProtein, (j->Position() > 0) ? 1 : 0, true);
                                outPeptideHeap.push_back(pepDecoy);
                                push_heap(outPeptideHeap.begin(), outPeptideHeap.end(),
                                                greater<TideIndexPeptide>());
                                ++decoysGenerated;
                        }
                }
        } else if (!allowDups) {
                for (set<string>::const_iterator i = setTargets.begin();
                                i != setTargets.end(); ++i) {
                        const string* setTarget = &*i;
                        const map<const string*, TargetInfo>::iterator targetLookup =
                                        targetInfo.find(setTarget);
                        const ProteinInfo& proteinInfo = (targetLookup->second.proteinInfo);
                        const int startLoc = targetLookup->second.start;
                        FLOAT_T pepMass = targetLookup->second.mass;
                        if (generateDecoy(*setTarget, targetToDecoy, &setTargets,
                                        &setDecoys, decoyType, allowDups, failedDecoyCnt,
                                        decoysGenerated, curProtein, proteinInfo, startLoc,
                                        pbProtein, pepMass, outPeptideHeap, outProteinSequences)) {

                                proteinWriter.Write(&pbProtein);
                        } else {
                                continue;
                        }
                }

        } else { // allow dups
                for (vector<pair<ProteinInfo, vector<PeptideInfo> > >::const_iterator i =
                                cleavedPeptideInfo.begin(); i != cleavedPeptideInfo.end();
                                ++i) {
                        for (vector<PeptideInfo>::const_iterator j = i->second.begin();
                                        j != i->second.end(); ++j) {
                                const string setTarget = j->Sequence();
                                const ProteinInfo& proteinInfo = i->first;
                                const int startLoc = j->Position();
                                FLOAT_T pepMass = calcPepMassTide(j->Sequence(), massType);
                                if (generateDecoy(setTarget, targetToDecoy, NULL, NULL,
                                                decoyType, allowDups, failedDecoyCnt, decoysGenerated,
                                                curProtein, proteinInfo, startLoc, pbProtein, pepMass,
                                                outPeptideHeap, outProteinSequences)) {
                                        proteinWriter.Write(&pbProtein);
                                } else {
                                        continue;
                                }
                        }
                }

        }
        if (invalidPepCnt > 0) {
                carp(CARP_INFO,
                                "Ignoring %d peptide sequences containing unrecognized characters",
                                invalidPepCnt);
        }
        if (failedDecoyCnt > 0) {
                carp(CARP_INFO,
                                "Failed to generate decoys for %d low complexity peptides",
                                failedDecoyCnt);
        }
        carp(CARP_INFO, "FASTA produced %d targets and %d decoys",
                        targetsGenerated, decoysGenerated);

        // Write to decoy fasta if necessary (if protein-reverse, we already wrote it)
        if (decoyFasta && decoyType != PROTEIN_REVERSE_DECOYS) {
                carp(CARP_INFO, "Writing decoy fasta...");
                // Iterate over all (protein, peptides from that protein)
                for (vector<pair<ProteinInfo, vector<PeptideInfo> > >::const_iterator i =
                                cleavedPeptideInfo.begin(); i != cleavedPeptideInfo.end();
                                ++i) {
                        string decoyProtein = *(i->first.sequence);
                        // Iterate over all peptides from the protein
                        for (vector<PeptideInfo>::const_iterator j = i->second.begin();
                                        j != i->second.end(); ++j) {
                                // In the protein sequence, replace the target peptide with its decoy
                                const string setTarget = j->Sequence();
                                const map<const string, const string*>::const_iterator decoyCheck =
                                                targetToDecoy.find(setTarget);
                                if (decoyCheck != targetToDecoy.end()) {
                                        decoyProtein.replace(j->Position(), j->Length(),
                                                        *(decoyCheck->second));
                                }
                        }
                        // Write out the final protein
                        (*decoyFasta) << ">" << decoyPrefix << i->first.name << endl
                                        << decoyProtein << endl;
                }
        }
}

void TideIndexApplication::writePeptidesAndAuxLocs(
                //vector<TideIndexPeptide>& peptideHeap,
                const string& peptidePbFile,

                const string& auxLocsPbFile, pb::Header& pbHeader) {
        // Check header
        if (pbHeader.source_size() != 1) {
                carp(CARP_FATAL, "pbHeader had a number of sources other than 1");
        }
        pb::Header_Source& headerSource = *(pbHeader.mutable_source(0));
        if (!headerSource.has_filename() || headerSource.has_filetype()) {
                carp(CARP_FATAL, "pbHeader source invalid");
        }

        string proteinsFile = headerSource.filename();
        vector<const pb::Protein*> proteins;
        pb::Header proteinsHeader;
        carp(CARP_INFO, "Reading proteins");
        if (!ReadRecordsToVector<pb::Protein>(&proteins, proteinsFile,
                        &proteinsHeader)) {
                carp(CARP_FATAL, "Error reading proteins from %s",
                                proteinsFile.c_str());
        } else if (proteinsHeader.file_type() != pb::Header::RAW_PROTEINS) {
                carp(CARP_FATAL, "Proteins file %s had invalid type",
                                proteinsFile.c_str());
        }
        // Clean up
        for (vector<const pb::Protein*>::iterator i = proteins.begin();
                        i != proteins.end(); ++i) {
                delete *i;
        }
        // The raw proteins file is read in. It's a valid source file;
        // remember it as such:
        headerSource.mutable_header()->CopyFrom(proteinsHeader);

        // Now check other desired settings
        if (!pbHeader.has_peptides_header()) {
                carp(CARP_FATAL, "!pbHeader->has_peptideHeapheader()");
        }
        const pb::Header_PeptidesHeader& settings = pbHeader.peptides_header();
        //if (!Peptide::SetMinMaxMassAndLength(settings)) {
        //  carp(CARP_FATAL, "Error setting min/max mass/length");
        if (!settings.has_enzyme() || settings.enzyme().empty()) {
                carp(CARP_FATAL, "Enzyme settings error");
        }

        pbHeader.set_file_type(pb::Header::PEPTIDES);
        pbHeader.mutable_peptides_header()->set_has_peaks(false);
        pbHeader.mutable_peptides_header()->set_decoys(
                        get_tide_decoy_type_parameter("decoy-format"));
        HeadedRecordWriter peptideWriter(peptidePbFile, pbHeader); // put header in outfile

        // Create the auxiliary locations header and writer
        pb::Header auxLocsHeader;
        auxLocsHeader.set_file_type(pb::Header::AUX_LOCATIONS);
        pb::Header_Source* auxLocsSource = auxLocsHeader.add_source();
        auxLocsSource->set_filename(peptidePbFile);
        auxLocsSource->mutable_header()->CopyFrom(pbHeader);
        HeadedRecordWriter auxLocWriter(auxLocsPbFile, auxLocsHeader);

        pb::Peptide pbPeptide;
        pb::AuxLocation pbAuxLoc;
        int auxLocIdx = -1;
        carp(CARP_DEBUG, "%d peptides in heap", outPeptideHeap.size());
        int count = 0;



        sort_heap(outPeptideHeap.begin(), outPeptideHeap.end(), greater<TideIndexPeptide>());


        while (!outPeptideHeap.empty()) {
                TideIndexPeptide curPeptide(outPeptideHeap.back());
                outPeptideHeap.pop_back();

                while (!outPeptideHeap.empty() && outPeptideHeap.back() == curPeptide){ //&& calcPepMassTide(peptideHeap.back().getSequence(), MONO) < 0) {
                      pb::Location* location = pbAuxLoc.add_location();
                      location->set_protein_id(outPeptideHeap.back().getProteinId());
                      location->set_pos(outPeptideHeap.back().getProteinPos());
                      outPeptideHeap.pop_back();
                }
                getPbPeptide(count, curPeptide, pbPeptide);
                // Not all peptides have aux locations associated with them. Check to seMONOe
                // if GetGroup added any locations to aux_location. If yes, only then
                // assign the corresponding array index to the peptide and write it out.
                if (pbAuxLoc.location_size() > 0) {
                        pbPeptide.set_aux_locations_index(++auxLocIdx);
                        auxLocWriter.Write(&pbAuxLoc);
                        pbAuxLoc.Clear();
                }

                // Write the peptide AFTER the aux_locations check, in case we added an
                // aux_locations_index to the peptide.
                peptideWriter.Write(&pbPeptide);

                if (++count % 100000 == 0) {
                        carp(CARP_INFO, "Wrote %d peptides", count);
                }
        }
}

FLOAT_T TideIndexApplication::calcPepMassTide(const string& sequence,
                MASS_TYPE_T massType) {
        FixPt mass;
        FixPt aaMass;
        if (massType == AVERAGE) {
                mass = MassConstants::fixp_avg_h2o;
                for (size_t i = 0; i < sequence.length(); ++i) {
                        if (i == 0) {
                                aaMass = MassConstants::fixp_nterm_avg_table[sequence[i]];
                        } else if (i == sequence.length() - 1) {
                                aaMass = MassConstants::fixp_cterm_avg_table[sequence[i]];
                        } else {
                                aaMass = MassConstants::fixp_avg_table[sequence[i]];
                        }
                        if (aaMass == 0) {
                                return -1;
                        }
                        mass += aaMass;
                }
        } else if (massType == MONO) {
                mass = MassConstants::fixp_mono_h2o;
                for (size_t i = 0; i < sequence.length(); ++i) {
                        if (i == 0) {
                                aaMass = MassConstants::fixp_nterm_mono_table[sequence[i]];
                        } else if (i == sequence.length() - 1) {
                                aaMass = MassConstants::fixp_cterm_mono_table[sequence[i]];
                        } else {
                                aaMass = MassConstants::fixp_mono_table[sequence[i]];
                        }
                        if (aaMass == 0) {
                                return -1;
                        }
                        mass += aaMass;
                }
        } else {
                carp(CARP_FATAL, "Invalid mass type");
        }
        return MassConstants::ToDouble(mass);
}

void TideIndexApplication::getPbProtein(int id, const string& name,
                const string& residues, pb::Protein& outPbProtein) {
        outPbProtein.Clear();
        outPbProtein.set_id(id);
        outPbProtein.set_name(name);
        outPbProtein.set_residues(residues);
}

void TideIndexApplication::getDecoyPbProtein(int id,
                const ProteinInfo& targetProteinInfo, string decoyPeptideSequence,
                int startLoc, pb::Protein& outPbProtein) {



        const string* proteinSequence = targetProteinInfo.sequence;

        const int pepLen = decoyPeptideSequence.length();

        // Add N term to decoySequence, if it exists
        if (startLoc > 0) {
                decoyPeptideSequence.insert(0, 1, proteinSequence->at(startLoc - 1));
        }
        // Add C term to decoySequence, if it exists
        size_t cTermLoc = startLoc + pepLen;
        decoyPeptideSequence.push_back(
                        (cTermLoc < proteinSequence->length()) ?
                                        proteinSequence->at(cTermLoc) : '-');
        // Append original target sequence, unless using protein level decoys
        if (get_tide_decoy_type_parameter("decoy-format")
                        != PROTEIN_REVERSE_DECOYS) {
                decoyPeptideSequence.append(
                                targetProteinInfo.sequence->substr(startLoc, pepLen));
        }

        getPbProtein(id, Params::GetString("decoy-prefix") + targetProteinInfo.name,
                        decoyPeptideSequence, outPbProtein);
        outPbProtein.set_target_pos(startLoc);
}

void TideIndexApplication::getPbPeptide(int id, const TideIndexPeptide& peptide,
                pb::Peptide& outPbPeptide) {
        outPbPeptide.Clear();
        outPbPeptide.set_id(id);
        outPbPeptide.set_mass(peptide.getMass());
        outPbPeptide.set_length(peptide.getLength());
        outPbPeptide.mutable_first_location()->set_protein_id(
                        peptide.getProteinId());
        outPbPeptide.mutable_first_location()->set_pos(peptide.getProteinPos());
        ;
        outPbPeptide.set_is_decoy(peptide.isDecoy());
}

void TideIndexApplication::addAuxLoc(int proteinId, int proteinPos,
                pb::AuxLocation& outAuxLoc) {
        pb::Location* location = outAuxLoc.add_location();
        location->set_protein_id(proteinId);
        location->set_pos(proteinPos);
}

void TideIndexApplication::processParams() {
        // Update mods-spec parameter for default cysteine mod
        string default_cysteine = "C+" + StringUtils::ToString(CYSTEINE_DEFAULT);
        string mods_spec = Params::GetString("mods-spec");
        if (mods_spec.find('C') == string::npos) {
                mods_spec =
                                mods_spec.empty() ?
                                                default_cysteine : default_cysteine + ',' + mods_spec;
                carp(CARP_DEBUG, "Using default cysteine mod '%s' ('%s')",
                                default_cysteine.c_str(), mods_spec.c_str());
        }
        Params::Set("mods-spec", mods_spec);
}

string getModifiedPeptideSeq(const pb::Peptide* peptide,
                const ProteinVec* proteins) {
        int mod_index;
        double mod_delta;
        stringstream mod_stream;
        const pb::Location& location = peptide->first_location();
        const pb::Protein* protein = proteins->at(location.protein_id());
        // Get peptide sequence without mods
        string pep_str = protein->residues().substr(location.pos(),
                        peptide->length());

        // Store all mod indices/deltas
        map<int, double> mod_map;
        set<int> mod_indices;
        for (int j = 0; j < peptide->modifications_size(); ++j) {
                //        var_mod_table.DecodeMod(ModCoder::Mod(peptide->modifications(j)), &mod_index, &mod_delta);
                MassConstants::DecodeMod(ModCoder::Mod(peptide->modifications(j)),
                                &mod_index, &mod_delta);
                mod_indices.insert(mod_index);
                mod_map[mod_index] = mod_delta;
        }
        for (set<int>::const_reverse_iterator j = mod_indices.rbegin();
                        j != mod_indices.rend(); ++j) {
                // Insert the modification string into the peptide sequence
                mod_stream << '[' << mod_map[*j] << ']';
                pep_str.insert(*j + 1, mod_stream.str());
                mod_stream.str("");
        }
        return pep_str;
}

bool TideIndexApplication::generateDecoy(const string& setTarget,
                std::map<const string, const string*>& targetToDecoy,
                set<string>* setTargets, set<string>* setDecoys, DECOY_TYPE_T decoyType,
                bool allowDups, unsigned int& failedDecoyCnt,
                unsigned int& decoysGenerated, int& curProtein,
                const ProteinInfo& proteinInfo, const int startLoc,
                pb::Protein& pbProtein, FLOAT_T pepMass,
                vector<TideIndexPeptide>& outPeptideHeap,
                vector<string*>& outProteinSequences) {

        const map<const string, const string*>::const_iterator decoyCheck =
                        targetToDecoy.find(setTarget);
        string* decoySequence = new string;
        if (decoyCheck != targetToDecoy.end()) {
                // Decoy already generated for this sequence
                *decoySequence = *(decoyCheck->second);
        } else {
                // Try to generate decoy
                if (allowDups) {
                        set<string> targets, decoys;
                        setTargets = &targets;
                        setDecoys = &decoys;
                }
                if (!GeneratePeptides::makeDecoy(setTarget, *setTargets, *setDecoys,
                                decoyType == PEPTIDE_SHUFFLE_DECOYS, *decoySequence)){
               /* if(!wmuMakeDecoy(setTarget, *setTargets, *setDecoys,
                                decoyType == PEPTIDE_SHUFFLE_DECOYS, *decoySequence, "NC")) {*/
                        carp(CARP_DETAILED_INFO, "Failed to generate decoy for sequence %s",
                                        setTarget.c_str());
                        ++failedDecoyCnt;
                        delete decoySequence;
                        return false;
                } else if (!allowDups) {
                        targetToDecoy[setTarget] =
                                        &*(setDecoys->insert(*decoySequence).first);
                } else {
                        targetToDecoy[setTarget] = decoySequence;
                }
        }

        outProteinSequences.push_back(decoySequence);

        // Write pb::Protein
        getDecoyPbProtein(++curProtein, proteinInfo, *decoySequence, startLoc,
                        pbProtein);
        // Add decoy to heap
        TideIndexPeptide pepDecoy(pepMass, setTarget.length(), decoySequence,
                        curProtein, (startLoc > 0) ? 1 : 0, true);
        outPeptideHeap.push_back(pepDecoy);
        push_heap(outPeptideHeap.begin(), outPeptideHeap.end(),
                        greater<TideIndexPeptide>());
        ++decoysGenerated;
        return true;
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */