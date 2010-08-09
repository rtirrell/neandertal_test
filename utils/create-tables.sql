CREATE TABLE neandertal_snps (
  rsid          VARCHAR(16) NOT NULL PRIMARY KEY,
  derived       VARCHAR(1)  NOT NULL,
  ancestral     VARCHAR(1)  NOT NULL,
  out_of_africa VARCHAR(1)  NOT NULL
);

