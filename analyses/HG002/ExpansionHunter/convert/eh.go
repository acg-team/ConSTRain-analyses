package main

import (
	"encoding/csv"
	"encoding/json"
	"flag"
	"fmt"
	"io"
	"log"
	"os"
)

var (
	colIdx = map[string]int{
		"chr":    0,
		"start":  1,
		"end":    2,
		"period": 3,
		"unit":   4,
	}
)

type jsonRecord struct {
	ID        string `json:"LocusID"`
	Structure string `json:"LocusStructure"`
	Region    string `json:"ReferenceRegion"`
	VarType   string `json:"VariantType"`
}

func main() {
	panelPtr := flag.String("i", "", "path to STR bed file to convert to ExpansionHunter json")

	flag.Parse()
	if *panelPtr == "" {
		log.Fatal("-i was not specified with no default")
	}

	f, err := os.Open(*panelPtr)
	if err != nil {
		log.Fatal(err)
	}
	defer f.Close()

	r := csv.NewReader(f)
	r.Comma = '\t'
	r.FieldsPerRecord = len(colIdx)

	// set capacity to slightly above known number of records in the bed file
	b := make([]jsonRecord, 0, 1.7e6)
	for {
		record, err := r.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			log.Fatal(err)
		}

		id := fmt.Sprintf("%s_%s", record[colIdx["chr"]], record[colIdx["start"]])
		structure := fmt.Sprintf("(%s)*", record[colIdx["unit"]])
		region := fmt.Sprintf("%s:%s-%s", record[colIdx["chr"]], record[colIdx["start"]], record[colIdx["end"]])

		jrecord := jsonRecord{ID: id, Structure: structure, Region: region, VarType: "Repeat"}
		b = append(b, jrecord)
	}

	jsonOut, err := json.MarshalIndent(
		b,
		"",
		"    ",
	)
	if err != nil {
		log.Fatal(err)
	}

	fmt.Println(string(jsonOut))
}
