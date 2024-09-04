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
	bufSize int = 1.7e6 // based on size of STR panel (1659608 lines)
)

type jsonRecord struct {
	ID        string `json:"LocusId"`
	Structure string `json:"LocusStructure"`
	Region    string `json:"ReferenceRegion"`
	VarType   string `json:"VariantType"`
}

func lineToJSON(record []string) jsonRecord {
	id := fmt.Sprintf("%s_%s", record[colIdx["chr"]], record[colIdx["start"]])
	structure := fmt.Sprintf("(%s)*", record[colIdx["unit"]])
	region := fmt.Sprintf("%s:%s-%s", record[colIdx["chr"]], record[colIdx["start"]], record[colIdx["end"]])

	return jsonRecord{ID: id, Structure: structure, Region: region, VarType: "Repeat"}
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

	// read lines and send them over `c`. `c` is buffered so we
	// should hopefully never block while reading
	c := make(chan []string, bufSize)
	go func() {
		// close channel once all lines are read
		defer close(c)

		for {
			record, err := r.Read()
			if err == io.EOF {
				break
			}
			if err != nil {
				log.Fatal(err)
			}

			c <- record
			// fmt.Println("Line read")
		}
	}()

	// as long as we're getting records over `c`, convert them to
	// `jsonRecord`s and add them to our buffer `b`
	b := make([]jsonRecord, 0, bufSize)
	for record := range c {
		b = append(b, lineToJSON(record))
		// fmt.Println("Line processed")
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
