Convert all KGML pathway files to an interactome. 

# Download KGML files 

## Find ID of organism
`curl -s  "http://rest.kegg.jp/list/organism" | grep "Homo sapiens"`

## Download pathway xml files
`curl "http://rest.kegg.jp/list/pathway/T01001" | cut -f 1 | while read A; do  curl -o "${A:5}.xml" "http://rest.kegg.jp/get/${A}/kgml" ; done`
