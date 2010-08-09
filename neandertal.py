import os
import sys
import subprocess

import genotype_tools
import mysql.connector

DB_INFO = {
  "host":     "marlowe",
  "user":     "gene210-user",
  "passwd":   "genomics",
  "buffered": True
}

# Valid population identifiers.
POPULATIONS = ["CEU", "YRI", "JPT", "CHB"]
POPULATION_CHOICES = ", ".join(POPULATIONS)

SNPS_PATH = "data/neandertal-SNPs.csv"
BASE_MAP = {
  "A": "T",
  "T": "A",
  "C": "G",
  "G": "C"
}

##
# Fetch information from the database for the user's value for a neandertal 
# SNP.
def get_snp_info(nsnp, flipped=False):
  # Score is the amount by which we increment the numerator and denominator of
  # the neandertal score.
  snp_info = {"score": (0, 0)}

  query = "SELECT * FROM diseases.neandertal_snps WHERE rsid = %s;"
  db.execute(query, (nsnp.rsid,))
  res = db.fetchone()

  # db.description is a list of tuples; the first element of each tuple is the 
  # field name.
  for i, attribute_info in enumerate(db.description):
    snp_info[attribute_info[0]] = res[i]

  snp_info["user_alleles"] = nsnp.genotype
  heterozygote = snp_info["out_of_africa"] + snp_info["ancestral"]

  # Homozygous homo neanderthalensis.
  if nsnp.genotype == 2 * snp_info["out_of_africa"]:
    snp_info["score"] = (2, 2)

  # Homozygous homo sapiens.
  elif nsnp.genotype == 2 * snp_info["ancestral"]:
    snp_info["score"] = (0, 2)

  # Otherwise, we're a heterozygote.
  elif sorted(nsnp.genotype) == sorted(heterozygote):
    snp_info["score"] = (1, 2)

  else:
    if not flipped:
      nsnp.genotype = "".join([BASE_MAP[a] for a in nsnp.genotype])
      get_snp_info(nsnp, True)
    else:
      # Hmmm... make no adjustment to the running numerator or denominator for 
      # this SNP, we couldn't match alleles. This ought never happen.
      print "Error occurred matching user SNP to",
      print "neandertal SNP %s!" % (nsnp.rsid,)

  return snp_info

##
# Point of entry is here!

if len(sys.argv) not in (2, 3):
  print "Usage: python neandertal.py path/to/user/genome/file.txt",
  print "[population identifier] (if population not given, you will be",
  print "prompted for its value."
  sys.exit(2)

population = None
if len(sys.argv) == 2:
  population = raw_input("Select a population (%s): " % POPULATION_CHOICES)
  population = population.strip().upper()

if len(sys.argv) == 3:
  population = sys.argv[2]

if population not in POPULATIONS:
  print "Population identifier '%s' is not valid." % (population, ),
  print "Valid population identifiers are one of (%s)." % (POPULATION_CHOICES,)
  sys.exit(2)

user_genome_path = sys.argv[1]
if not os.path.exists(user_genome_path):
  print "No file found at '%s'." % (user_genome_path)
  sys.exit(2)

print "Loading genotype file '%s'... " % (user_genome_path,),
sys.stdout.flush()
user_snps = genotype_tools.FileUtils.read_genotype_file(user_genome_path)
print "done!"

##
# Calculation begins at this point.
results = {}

database = mysql.connector.Connect(**DB_INFO)
db = database.cursor()
db.execute("SELECT rsid FROM diseases.neandertal_snps;")

for nrsid in [row[0] for row in db.fetchall()]:
  user_nsnp = user_snps.get(nrsid, None)
  
  # We could not get this SNP directly: impute.
  if user_nsnp is None:
    user_nsnp = genotype_tools.impute_rsid_simple(user_snps, nrsid, population)
    print "Imputed %s -> %s." % (user_nsnp.nearest_SNP, user_nsnp.rsid,)

  snp_info = get_snp_info(user_nsnp)
  # If we imputed, user_nsnp.rsid is the imputed rsid, and nearest_SNP is the
  # rsid of the SNP we were imputing from (imputing for).
  if user_nsnp.nearest_SNP is not None:
    snp_info["imputed_from"] = user_nsnp.nearest_SNP
  results[user_nsnp.rsid] = snp_info

nnumerator   = sum(a["score"][0] for a in results.values())
ndenominator = sum(a["score"][1] for a in results.values())

# Name the output HTML file based on the input genome file"s name.
filename = os.path.splitext(os.path.split(user_genome_path)[1])[0]
out_file_name = "%s.html" % (filename,)
out_file = open(out_file_name, "w")

out_file.write("""
<html>
  <head>
    <title>
      Are You a Caveman?
    </title>
  </head>
  <body>
  <h3 style="text-align:left;">
    Caveman Index:
  </h3>
  Genome file: %s<br>
  Population given was: %s<br>
  Alleles shared with a caveman: %d / %d (%0.4f%%)<br><br>
""" % (filename, population, nnumerator, ndenominator, 
       100 * nnumerator / float(ndenominator)))

#chart_url = build_chart_url(el_running_odds)
#out_file.write("<img src='%s'></img><br><br>" % (chart_url,))

out_file.write("""
  <table border="1">
    <tr>
      <th>rsID</th>
      <th>Caveman Allele</th>
      <th>Ancestral Allele</th>
      <th>User Alleles</th>
      <th>Imputed From</th>
    </tr>
""")

# We just reserve the right to sort later, by changing the key function to
# something meaningful (this code was extensively copied and pasted).
for (rsid, snp_info) in sorted(
  results.items(), key=lambda v: v, reverse=True
):
  out_file.write("""
    <tr>
      <td>%s</td>
      <td>%s</td>
      <td>%s</td>
      <td>%s</td>
  """ % (rsid, snp_info["out_of_africa"], 
         snp_info["ancestral"], "/".join(snp_info["user_alleles"])))

  if "imputed_from" in snp_info:
    out_file.write("<td>%s</td>" % (snp_info["imputed_from"]))

  out_file.write("</td>")

out_file.write("</table></body></html>")
out_file.close()
subprocess.Popen(("open", out_file_name)).wait()
