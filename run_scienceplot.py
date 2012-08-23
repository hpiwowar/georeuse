import pprint
import os
import time
import re
import nose
from nose.tools import assert_equals
from utils.testing import slow, online
from utils.tidbits import flatten_unique
from collections import defaultdict
import EUtils
from EUtils import HistoryClient, ThinClient
import dataset
import datasources
from datasources import pubmedcentral
from datasources import pubmed
from datasources import affiliation
from datasources import geo
from datasources import geo_reuse
from datasources import oaexcerpt
from datasources import urlopener
from utils.cache import TimedCache
import pickle
import csv

EMAIL_CONTACT = "hpiwowar@gmail.com"
VERBOSE = False

#base_query = """(GEO[text] OR omnibus[text]) NOT "pmc gds"[filter]"""
#base_query = """(GEO[text] OR omnibus[text]) NOT "pmc gds"[filter] AND ("1900"[PubDate] : "2009"[PubDate])"""
base_query_reuse = """("1900"[PubDate] : "2012/06/30"[PubDate]) NOT "pmc gds"[filter]"""
base_query_submit = """("1900"[PubDate] : "2012/06/30"[PubDate]) AND "pmc gds"[filter]"""

#url_for_gse = """http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term=GSE[ETYP]&retmax=10000&usehistory=n"""
#url_for_gds = """http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term=GDS[ETYP]&retmax=10000&usehistory=n"""

@slow
def test_get_accession_in_pmc_fulltext():
    response = get_accession_in_pmc_fulltext("GSE", "200008514", base_query_reuse)
    assert_equals(response, ['2785812', '2620272'])

    response = get_accession_in_pmc_fulltext("GSE", "200008478", base_query_reuse)
    assert_equals(response, ['2223695'])

    response = get_accession_in_pmc_fulltext("GSE", "200007778", base_query_reuse)
    assert_equals(response, [])

    response = get_accession_in_pmc_fulltext("GDS", "2643", base_query_reuse)
    assert_equals(response, ['2715883', '2781753'])
    
@slow    
def test_get_accession_variants():
    response = get_accession_variants("GSE", "200008514")
    assert_equals(response, ['GSE200008514', '"GSE 200008514"', 'GSE8514', '"GSE 8514"'])
    
def get_accession_variants(id_type, id):
    variants = []
    variants.append(id_type + id)
    variants.append('"' + id_type + ' ' + id  + '"')
    id_stripped = geo.get_stripped_accession(id)
    if (id_stripped != id):
        variants.append(id_type + id_stripped)
        variants.append('"' + id_type + ' ' + id_stripped + '"')
    return(variants)
        
def get_accession_in_pmc_fulltext(id_type, id, pmc_query):
    pmc_ids = []
    accession_variants = get_accession_variants(id_type, id)
    for variant in accession_variants:
        query = variant + "[text] AND " + pmc_query
        try:
            pmc_ids += pubmedcentral.search(query)        
        except Exception, e:
            print "error for some reason!", e
            print "skipping query= ", query
            print 1/0
            pass
    return(pmc_ids)
    
@slow    
def test_get_authors_and_submittors_from_accession():
    response = get_authors_and_submittors_from_accession([u'17911395'], [u'Ye P', u'Rainey W'])
    assert_equals(response, ['Mariniello', 'Shibata', 'Ye', 'Rainey', 'Mantero'])

def get_author_last_names(pmids):
    last_names_list = []
    for pmid in pmids:
        authors = pubmed.authors([pmid])
        authors_list = authors[0].split(";")
        last_names = [author.split(" ")[0] for author in authors_list]
        last_names_list += last_names
    return(last_names_list)
    
def get_authors_and_submittors_from_accession(pmids, contributors=None):
    last_names_list = get_author_last_names(pmids)
    if contributors:
        contributor_last_names = [author.split(" ")[0] for author in contributors]
        last_names_list += contributor_last_names
    last_names_set = set(last_names_list)
    return(list(last_names_set))
    
    
@slow
def test_get_author_intersect_submit_reuse():
    response = get_author_intersect_submit_reuse(['Pessi', 'Ahrens', 'Rehrauer', 'Lindemann', 'Hauser', 'Fischer', 'Hennecke', 'Lindemann', 'Moser', 'Pessi', 'Hauser', 'Friberg', 'Hennecke', 'Fischer', u'Pessi', u'Ahrens', u'Rehrauer', u'Lindemann', u'Hauser', u'Fischer', u'Hennecke'], ['Hacker', 'Sohlenkamp', 'Pessi', 'Aktas', 'Geiger', 'Narberhaus'])
    assert_equals(response, ['Pessi'])
    
def get_author_intersect_submit_reuse(submit_authors, reuse_authors):
    intersect = multi_intersection([submit_authors, reuse_authors])
    return(intersect)

meshes = """"Algorithms"[mesh]
"Databases, Genetic"[mesh]
"Gene Expression Profiling/methods"[mesh]
"Computational Biology/methods"[mesh]
"Oligonucleotide Array Sequence Analysis/methods"[mesh]
"Genomics/methods"[mesh]
"Reproducibility of Results"[mesh]
"Software"[mesh]
"Computer Simulation"[mesh]
"Internet"[mesh]
"Data Interpretation, Statistical"[mesh]""".split("\n")

metaquery = "(meta-analysis [pt] OR meta-analysis [tw] OR metanalysis [tw]) OR meta-analysis [mh])"

words = """submitted
deposited
user*
public
accessed
downloaded
published""".split("\n")


    
def test_get_all_attributes():
    test_dict = {('GSE8514', '2785812'): ('GSE8514', [u'17911395'], '2785812', [u'19917117']), ('GSE8514', '2620272'): ('GSE8514', [u'17911395'], '2620272', [u'19014681']), ('GSE8478', '2223695'): ('GSE8478', [u'17977147', u'17951393'], '2223695', [u'17993534'])}
    response = get_all_attributes(test_dict)
    assert_equals(response, "hi heather")
   
def get_all_GSE_GDS_excerpts(pmcids):
    pattern = "GSE|GDS"
    excerpts = [oaexcerpt.get_oa_excerpt(pmid, pattern) for pmcid in pmcids]
    return(excerpts)
    
def get_all_attributes(id_dict):
    all_accession_key = [vals[0] for vals in id_dict.values()]
    all_submit_pmids = flatten_unique([vals[3] for vals in id_dict.values()])
    all_reuse_pmcids = flatten_unique([vals[4] for vals in id_dict.values()])
    all_reuse_pmids = flatten_unique([vals[5] for vals in id_dict.values()])

    reuse_affiliation = affiliation.institution(all_reuse_pmids)
    journal         = pubmed.journal(all_reuse_pmids)
    year            = pubmed.year_published(all_reuse_pmids)
    date_published  = pubmed.date_published(all_reuse_pmids)
    medline_status  = pubmed.medline_status(all_reuse_pmids)
    is_geo_reuse    = geo_reuse.is_geo_reuse(all_reuse_pmids)
    reuse_is_oa     = pubmed.is_open_access(all_reuse_pmids)
    metaanal        = pubmed._get_flags_for_pattern(all_reuse_pmids, metaquery)
    
    oa_excerpts     = oaexcerpt.get_oa_excerpts(all_reuse_pmcids, "(GSE.\d|GDS.\d|omnibus)", flags = re.IGNORECASE|re.MULTILINE)

    biolink_filter   = pubmedcentral.get_flags_for_pattern(all_reuse_pmcids, '(geo OR omnibus)  AND microarray  AND "gene expression" AND accession NOT (databases OR user OR users  OR (public AND accessed) OR (downloaded AND published))')
    basic_reuse_filter = pubmedcentral.get_flags_for_pattern(all_reuse_pmcids, '"gene expression omnibus" AND (submitted OR deposited)')
    creation_filter    = pubmedcentral.get_flags_for_pattern(all_reuse_pmcids, '("gene expression" [text] AND "microarray" [text] AND "cell" [text] AND "rna" [text]) AND ("rneasy" [text] OR "trizol" [text] OR "real-time pcr" [text]) NOT ("tissue microarray*" [text] OR "cpg island*" [text])')

    has_mesh = {}
    for term in meshes:
        has_mesh[term] = pubmed.filter_pmids(all_reuse_pmids, term)
    mesh_filters = [";".join([term for term in has_mesh if pmid in has_mesh[term]]) for pmid in all_reuse_pmids]
    
    has_word = {}
    for word in words:
        has_word[word] = pubmedcentral.filter_pmcids(all_reuse_pmcids, word)
    word_filters = [";".join([word for word in has_word if pmcid in has_word[word]]) for pmcid in all_reuse_pmcids]
    print word_filters
    
    reuse_pmid_dict = defaultdict(tuple, zip(all_reuse_pmids, zip(reuse_affiliation, journal, year, date_published, medline_status, is_geo_reuse, reuse_is_oa, metaanal, mesh_filters)))
    reuse_pmcid_dict = defaultdict(tuple, zip(all_reuse_pmcids, zip(biolink_filter, basic_reuse_filter, creation_filter, oa_excerpts, word_filters)))

    full_dict = {}
    for vals in id_dict.values():
        id = vals[0]
        reuse_pmcid = vals[4]
        reuse_pmid = vals[5][0] if vals[5] else ""
        full_dict[id+reuse_pmcid] = vals + ("|",) + reuse_pmcid_dict[reuse_pmcid] + ("|",) + reuse_pmid_dict[reuse_pmid]
        print full_dict[id+reuse_pmcid]
    return(full_dict)
    
def test_get_response_dict():
    response_dict = get_response_dict("GSE", ["200008514", "200008478"], base_query_reuse)
    assert_equals(response_dict.items()[0:2], [(('GSE8514', '2785812'), ('GSE8514', ['8514'], ['2860'], [u'17911395'], '2785812', [u'19917117'], ['Mariniello', 'Shibata', 'Ye', 'Rainey', 'Mantero'], ['Kroll', 'Barkema', 'Carlon'], [], ['Medical College of Georgia'], u'Jul 19, 2007')), (('GSE8514', '2620272'), ('GSE8514', ['8514'], ['2860'], [u'17911395'], '2620272', [u'19014681'], ['Mariniello', 'Shibata', 'Ye', 'Rainey', 'Mantero'], ['Tozeren', 'Gormley'], [], ['Medical College of Georgia'], u'Jul 19, 2007'))])
    
def get_response_dict(id_type, ids, pmc_query):
    response_dict = {}
    num_ids = len(ids)
    id_counter = 0
    geo_instance = geo.GEO()
    for accession in ids:
        print accession,
        id_counter += 1
        reuse_pmcids = get_accession_in_pmc_fulltext(id_type, accession, pmc_query)        
        if ((not reuse_pmcids) or ("<ERROR>" in "".join(reuse_pmcids))):
            print " Nope"
            continue
        stripped_accession = geo.get_stripped_accession(accession)  
        if id_type=="GSE":
            gse_accessions = [stripped_accession]
            gds_accessions = [geo.get_stripped_accession(acc) for acc in geo.get_gds_from_gse("GSE"+stripped_accession)]
        else:
            gds_accessions = [stripped_accession]  
            gse_accessions = [geo.get_stripped_accession(acc) for acc in geo.get_gse_from_gds("GDS"+stripped_accession)]
            
        try:       
            submit_pmids = geo_instance.pmids("GSE"+gse_accessions[0])
        except Exception:
            continue
            
        print id_counter, "of", num_ids, ":", stripped_accession, "--", (submit_pmids), "; ", len(reuse_pmcids)
        for reuse_pmcid in reuse_pmcids:
            reuse_pmids_for_pmcid = pubmedcentral.pmcids_to_pmids(reuse_pmcid)
            dict_key = (id_type+stripped_accession, reuse_pmcid)
            this_submit_contributors = flatten_unique([geo_instance.contributors("GSE"+gse_accession) for gse_accession in gse_accessions])
            this_submit_authors = get_authors_and_submittors_from_accession(submit_pmids, this_submit_contributors)
            this_reuse_authors   = get_authors_and_submittors_from_accession(reuse_pmids_for_pmcid)
            intersect = get_author_intersect_submit_reuse(this_submit_authors, this_reuse_authors) 
            submit_affiliation = affiliation.institution(submit_pmids)
            release_date = geo_instance.release_date("GSE"+gse_accession)     
            response_dict[dict_key] = (id_type+stripped_accession, gse_accessions, gds_accessions, submit_pmids, reuse_pmcid, reuse_pmids_for_pmcid, this_submit_authors, this_reuse_authors, intersect, submit_affiliation, release_date)
            #print response_dict[dict_key]
    return(response_dict)

def test_get_dict_submit_reuse_of_accession_in_pmc_fulltext():
    full_dict = get_dict_submit_reuse_of_accession_in_pmc_fulltext("GSE", ["200008514", "200008478"], base_query_reuse)
    assert_equals(full_dict.items()[0:2], [('GSE85142620272', ('GSE8514', ['8514'], ['2860'], [u'17911395'], '2620272', [u'19014681'], ['Mariniello', 'Shibata', 'Ye', 'Rainey', 'Mantero'], ['Tozeren', 'Gormley'], [], ['Medical College of Georgia'], u'Jul 19, 2007', '|', '1', '0', '0', 'issue Phenotype Data Tissue No. of Samples Gene Expression Omnibus/Array Express Accn. # Adipose 10 GSE3526 Adrenal 20 GSE3526, GSE8514, GSE2316 Brain 89 GSE3526, GSE7621, GSE7307, GSE2361, E_AFMX-11, E-TA|20, Colon 10 E-TABM-176, GSE8671, GSE9254, GSE9452 Epidermal 25 GSE1133, GSE2361, GSE3419, GSE3526, GSE7307 Heart 38 E_AFMX-11, E-MIMR-27, GSE1133, GSE2240, GSE2361, GSE3526, GSE3585, GSE7307 Kidney 10 E_A|SE2004, GSE2361, GSE3526, GSE7392 Liver 10 E_AFMX-11, GSE2004, GSE3526, GSE6764 Lung 26 E-MEXP-231, GSE10072, GSE1133, GSE2361, GSE3526 Mammary 15 E-TABM-66, GSE2361, GSE3526, GSE7307, GSE7904 Muscle 64 GS|GSE3526, GSE5110, GSE6798, GSE7307, GSE9103, Ovary 10 GSE2361, GSE3526, GSE6008, GSE7307 Pancreas 6 GSE1133, GSE2361, GSE7307 Peripheral blood 12 GSE7462, GSE8608, GSE8668, GSE8762, GSE9692 Small intestine|een 12 GSE2004, GSE2361, GSE3526, GSE7307 Stomach 10 GSE2361, GSE3526, GSE7307 Testis 38 E_AFMX-11, GSE1133, GSE2361, GSE3218, GSE3526, GSE7307, GSE7808 Thymus 5 GSE1133, GSE2361, GSE7307 Infectious Diseas|. of Samples Gene Expression Omnibus/Array Express Accn. # Hepititis C 147 GSE11190, GSE7123 HIV 41 GSE6740, GSE9927 Influenza A 28 GSE6269 Malaria 15 GSE5418 Table 2 Adjusted Rand Index compares observed ', '|', 'Drexel University', 'BMC Bioinformatics', '2008', '2008 Nov 17', 'indexed for MEDLINE', '1', '1', '0', '"Algorithms"[mesh]')), ('GSE85142785812', ('GSE8514', ['8514'], ['2860'], [u'17911395'], '2785812', [u'19917117'], ['Mariniello', 'Shibata', 'Ye', 'Rainey', 'Mantero'], ['Kroll', 'Barkema', 'Carlon'], [], ['Medical College of Georgia'], u'Jul 19, 2007', '|', '1', '0', '0', '', '|', 'Katholieke Universiteit Leuven', 'Algorithms Mol Biol', '2009', '2009 Nov 16', 'PubMed', '1', '1', '0', '"Computational Biology/methods"[mesh]'))])
        
def get_dict_submit_reuse_of_accession_in_pmc_fulltext(id_type, ids, pmc_query):
    response_dict = get_response_dict(id_type, ids, pmc_query)
    full_dict = get_all_attributes(response_dict)
    return(full_dict)
    
@slow    
def test_authors_in_common_from_pmids():
    pmids = ["20349403", "18998887", "18767901"]
    response = authors_in_common_from_pmids(pmids)
    assert_equals(response, ["Piwowar"])

def multi_intersection(xs):
    inter = reduce(set.intersection, [set(x) for x in xs])
    return list(inter)
    
def authors_in_common_from_pmids(pmids):
    last_names_list = []
    for pmid in pmids:
        authors = pubmed.authors([pmid])
        authors_list = authors[0].split(";")
        last_names = [author.split(" ")[0] for author in authors_list]
        last_names_list.append(last_names)
    authors_intersection = multi_intersection(last_names_list)
    return(authors_intersection)
    
def get_from_query_gds_in_pmc_fulltext_dict(id_type, pmc_query, geo_year=None):  
    ids = geo.get_ids_by_year(id_type, geo_year)
    ids.sort()
    response_dict = get_dict_submit_reuse_of_accession_in_pmc_fulltext(id_type, ids, pmc_query)  
    return(response_dict)


def print_accession_pmcids(prefix, accession_dict):
    for accession in accession_dict:
        pmcid_list = accession_dict[accession]
        if pmcid_list:
            for pmcid in pmcid_list:
                print "%s\t%s%s\t%s" %(prefix, prefix, accession, pmcid)
        else:
            print "%s\t%s%s\t%s" %(prefix, prefix, accession, "")

#print_accession_pmcids("GDS", gds_dict)
#print_accession_pmcids("GSE", gse_dict)

@slow
def test_estimate_pmc_coverage():
    response = estimate_pmc_coverage("arrayexpress[title]")
    assert_equals(response, (6, 14, 0.42857142857142855))

    # changes with time
    response = estimate_pmc_coverage('"gene expression profiling"[mesh]')
    #assert_equals(response, (11038, 47384, 0.23294783049130507))

@slow
def test_estimate_pmc_coverage_given_years():
    # changes with time
    response = estimate_pmc_coverage('"gene expression profiling"[mesh]', "2007", "2009")
    #assert_equals(response, (6311, 21569, 0.29259585516250175))
    
def estimate_pmc_coverage(query, start_year="1800", end_year="3000"):
    pubmed_query = query + ' AND ("' + start_year + '"[pdat] : "' + end_year + '"[pdat])'
    pubmed_ids = pubmed.search(pubmed_query)
    num_pubmed = len(pubmed_ids)

    pmc_query = query + ' AND ("' + start_year + '"[PubDate] : "' + end_year + '"[PubDate])'
    pmc_ids = pubmedcentral.search(pmc_query)
    num_pmc = len(pmc_ids)
    
    ratio = num_pmc / (num_pubmed + 0.0)
    
    return(num_pmc, num_pubmed, ratio)


fields = "accession	gse	gds	submit_pmids	reuse_pmcid	reuse_pmids_for_pmc	this_submit_authors	this_reuse_authors	intersect	submit_affiliation	release_date	sep1	bioloink_filter	basic_reuse_filter	creation_filter	oa_excerpts	word_filters	sep2	reuse_affiliation	journal	year	date_published	medline_status	is_geo_reuse	reuse_is_oa	metaanal	mesh_filters".split("\t")

def get_excerpts():
    lines = open("scienceplot_new/results/reuse_pmcids.txt").readlines()
    accessions = [line.split("\t")[0].strip() for line in lines]
    reuse_pmcids = [line.split("\t")[1].strip() for line in lines]
    cached_article_dir = "scienceplot_new/articles"
    oa_excerpts     = oaexcerpt.get_oa_excerpts(reuse_pmcids, "(GSE.\d|GDS.\d|omnibus|download|publicly)", 200, 200, re.IGNORECASE|re.MULTILINE, cached_article_dir)
    excerpts_file = open("scienceplot_new/results/excerpts.txt", "w")
    for (accession, pmcid, excerpt) in zip(accessions, reuse_pmcids, oa_excerpts):
        lookup_number = accession[3:]  # remove the prefix
        excerpt_tagged = re.sub(lookup_number, lookup_number + "{{tag}}", excerpt)        
        excerpts_file.write(accession + "\t" + pmcid + "\t"+ str(excerpt) + "\t" + excerpt_tagged + "\r\n")
    excerpts_file.close()

def run_stats(geo_years, id_types, base_query):
    response_dict = {}
    for geo_year in geo_years:
        print "\n\n****\n", geo_year
        for id_type in id_types:
            response_dict[(id_type, geo_year)] = get_from_query_gds_in_pmc_fulltext_dict(id_type, base_query, geo_year)
            pkl_file = open("scienceplot_new/results/" + id_type + "_dict" + geo_year + "b.pkl", "wb")
            pickle.dump(response_dict[(id_type, geo_year)], pkl_file)
            pkl_file.close()
    
            fh = open("scienceplot_new/results/" + id_type + geo_year + "b.csv", "w")
            header = ",".join(fields)
            fh.write(header + "\n")
            dataset.csv_write_to_file(fh, response_dict[(id_type, geo_year)].values())
            fh.close()

    return(response_dict)

geo_years = [str(year) for year in range(2000,2012)]
id_types = ["GSE", "GDS"]

fh = open("scienceplot_new/results/pubmed_pmc_ratios.csv", "w")
ratio_fields = ["year", "num_pmc", "num_pubmed", "pmc_pmid_ratio"]
writer = csv.DictWriter(fh, ratio_fields)
writer.writerow(dict((fn,fn) for fn in writer.fieldnames))
for year in geo_years:
    (num_pmc, num_pubmed, ratio) = estimate_pmc_coverage('"gene expression profiling"[mesh]', year, year)
    row_dict = dict(year=year, num_pmc=num_pmc, num_pubmed=num_pubmed, pmc_pmid_ratio=round(100*ratio))
    print row_dict.values()
    writer.writerow(row_dict)
fh.close()

updated_base_query_all = """("1901"[PubDate] : "2011/12/31"[PubDate])"""
updated_base_query_reuse = """("1901"[PubDate] : "2011/12/31"[PubDate]) NOT "pmc gds"[filter]"""
updated_base_query_create = """("1901"[PubDate] : "2011/12/31"[PubDate]) AND "pmc gds"[filter]"""

if True:
    fh = open("scienceplot_new/results/pubmed_gse_count.csv", "w")
    num_id_fields = ["year", "num_gse_ids"]
    writer = csv.DictWriter(fh, num_id_fields)
    writer.writerow(dict((fn,fn) for fn in writer.fieldnames))
    total_count = 0
    for id_type in ["GSE"]:
        id_count = 0
        for year in geo_years:
            ids = geo.get_ids_by_year(id_type, year)
            row_dict = dict(year=year, num_gse_ids=len(ids))
            print row_dict.values()
            writer.writerow(row_dict)
            id_count += len(ids)
        print id_type, id_count
        total_count += id_count
    print total_count  
    fh.close()


year = "2009"
#submission_years = [str(year) for year in range(2007,2008)]
submission_years = [year]
all_dict = run_stats(submission_years, id_types, updated_base_query_all)

if True:
    all_dict_dicts = []
    for cohort in all_dict.keys():
        for acc in all_dict[cohort]:
            record = all_dict[cohort][acc]
            record_dict = dict(zip(fields, record))
            lookup = record_dict["accession"]
            lookup_number = lookup[3:]
            excerpt = record_dict["oa_excerpts"]
            excerpt_tagged = re.sub(lookup_number, lookup_number + "{{tag}}", excerpt)
            record_dict["tagged_oa_excerpts"] = excerpt_tagged
            record_dict["cohort"] = cohort
            all_dict_dicts.append(record_dict)
    new_fields = fields
    new_fields.append("tagged_oa_excerpts")
    new_fields.append("cohort")
    
    fh = open("scienceplot_new/results/records_" + year + ".csv", "w")
    writer = csv.DictWriter(fh, new_fields)
    for row in all_dict_dicts:
        writer.writerow(row)
    fh.close()
    
        

