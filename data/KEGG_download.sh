# Find ID of organism
curl -s  "http://rest.kegg.jp/list/organism" | grep "Arabidopsis thaliana"

# Download pathway xml files
curl "http://rest.kegg.jp/list/pathway/T01001" | cut -f 1 | while read A; do  curl -o "human_KGML/${A:5}.xml" "http://rest.kegg.jp/get/${A}/kgml" ; done

# Download compound ID mapping
curl "http://rest.kegg.jp/list/compound" > humanKGML/KEGG_compound_ids.txt