# This is a script for awk. To run: "awk -f pbrename.awk <input_dict.file> <replace.file>

NR == FNR {rep[$1] = $2 next} 
{
    for (key in rep) 
        gsub(key, rep[key]) 
    print
}
