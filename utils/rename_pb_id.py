#!/usr/bin/env python
# %%
from dttimeframe.timeFrame import tag
import argparse
# %%
parser = argparse.ArgumentParser(description="rename PB transcript id")
parser.add_argument("input", help="input filename", type=str)
parser.add_argument("output", help="output filename", type=str)
parser.add_argument("-c","--combine", help="output combination", type=str, default="")
args = parser.parse_args()
# %%
Tool = tag()
Tool.start()
comnbine_set = set()
from_str = ""
to_str = ""
target_list = ["gene_name","ref_gene_id","cmp_ref_gene"]
with open(args.input,"r") as target, open(args.output,"w") as output:
    for line in target:
        line_list = line.split("\t")
        type_str = line_list[2]
        if type_str == "transcript":
            attribute_str = line_list[8]
            attribute_list = attribute_str.split("; ")
            attribute_split_list = [n.replace("\"","").split(" ") for n in attribute_list]
            attribute_dict = {}
            for attribute_split in attribute_split_list:
                if len(attribute_split) > 2:
                    Tool.timeStamp("Error: \"{}\"".format(" ".join(attribute_split)))
                attribute_value = " ".join(attribute_split[1:])
                attribute_dict[attribute_split[0]] = attribute_value
            comnbine_set.update({" ".join(sorted(list(attribute_dict.keys())))})
            name_list = [attribute_dict.get(n,"") for n in target_list if attribute_dict.get(n,"") != ""]
            if len(name_list) > 0:
                from_str = "transcript_id \"{}\"".format(attribute_dict["transcript_id"])
                name_str = name_list[0]
                to_str = "transcript_id \"{}_{}\"".format(attribute_dict["transcript_id"],name_str)
                output.write(line.replace(from_str,to_str))
            else:
                from_str = ""
                to_str = ""
                output.write(line)
        elif type_str == "exon":
            if from_str == "":
                output.write(line)
            else:
                output.write(line.replace(from_str,to_str))
        else:
            output.write(line)
            Tool.timeStamp(F"Error: {line}")
# %%
if args.combine != "":
    with open(args.combine,"w") as combine_out:
        combine_out.write("\n".join(list(comnbine_set))+"\n")
# %%
Tool.stop()
