#!/bin/env ruby
# == NAME
# samplesheet_from_folder.rb
#
# == USAGE
# ./this_script.rb [ -h | --help ]
#[ -i | --infile ] |[ -o | --outfile ] | 
# == DESCRIPTION
# A script to produce a basic sample sheet for exome processing from a folder of fastQ files
#
# == OPTIONS
# -h,--help Show help
# -i,--infile=INFILE input file
# -o,--outfile=OUTFILE : output file

#
# == EXPERT OPTIONS
#
# == AUTHOR
#  Marc Hoeppner, mphoeppner@gmail.com

require 'optparse'
require 'ostruct'


### Define modules and classes here


### Get the script arguments and open relevant files
options = OpenStruct.new()
options.pacbio = false
opts = OptionParser.new()
opts.banner = "Reads Fastq files from a folder and writes a sample sheet to STDOUT"
opts.separator ""
opts.on("-f","--folder", "=FOLDER","Folder to scan") {|argument| options.folder = argument }
opts.on("-p", "--[no-]pacbio [FLAG]", TrueClass, "Pacbio data") {|argument| options.pacbio = argument.nil? ? true : argument }
opts.on("-c","--centre", "=CENTRE","Name of sequencing centre") {|argument| options.centre = argument }
opts.on("-s","--platform","=PLATFORM", "Name of the sequencing platform") {|argument| options.platform = argument }
opts.on("-l","--lookup", "=LOOKUP", "Lookup file with lib id <> other id") {|argument| options.lookup = argument }
opts.on("-h","--help","Display the usage information") {
 puts opts
 exit
}

opts.parse! 

abort "Folder not found (#{options.folder})" unless File.directory?(options.folder)

date = Time.now.strftime("%Y-%m-%d")
options.centre ? center = options.centre : center = "IKMB"
options.platform ? platform = options.platform : platform = "NovaSeq6000"

lookup = {}
if options.lookup
	IO.readlines(options.lookup).each do |line|
		key,value = line.strip.split("\t")
		lookup[key] = value
	end
end

if options.pacbio

	fastq_files = Dir["#{options.folder}/*.fastq.gz"]

	puts "patient;sample;R1"

	fastq_files.each do |file|

		name = file.gsub(".fastq.gz", "").split("/")[-1]
		path = File.expand_path(file)
		puts "#{name};#{name};#{path}"

	end

else 

	# 220600000039-A04_22Jun39-A04-L1_S1_L001_R1_001.fastq.gz

	fastq_files = Dir["#{options.folder}/*_R*.fastq.gz"].sort

	groups = fastq_files.group_by{|f| f.split("/")[-1].split(/_R[1,2]/)[0] }

	puts "patient;sample;library;rgid;rgpu;R1;R2"

	groups.each do |group, files|

            left,right = files.sort.collect{|f| File.absolute_path(f)}

            abort "Missing one member of the pair for #{group}" unless left && right

            library = group.split("_")[1]
            sample = library

            if lookup.has_key?(library)
                sample = "#{lookup[library]}"
            end

            e = `zcat #{left} | head -n1 `
            e.gsub!("@", "")
            header = e

            instrument,run_id,flowcell_id,lane,tile,x,y = header.split(" ")[0].split(":")

            index = header.split(" ")[-1].split(":")[-1]
            readgroup = flowcell_id + "." + lane + "." + library 

            pgu = flowcell_id + "." + lane + "." + index

            puts "#{sample};#{sample};#{library};#{readgroup};#{pgu};#{left};#{right}"

	end

end


