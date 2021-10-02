# Kc-Nanoindentation-Analysis-software
Youtube video tutorial: 

Main purpose: to extract the work from load-control nanoindentation load-unload curves from , fit a power-law to the curves in order to calculate the fracture toughness without imaging the indents.

File format: .txt file. See file_format_example_01.jpg for an example of the file format. Let's say you have 10 indents load-displacement curves you want to analyze. The .txt file should have 20 columns in total. Each pair of columns = 1 indent measurement displacement and load, respectively. IMPORTANT: The software does not work if your data has an uneven number of rows and it also doesn't work if you have a blank column somewhere. So use Microsoft Excel (or your preferred excel-like software) to delete blank columns and to trim the bottom of your data to make sure all the columns have the same # of rows. See file_format_example_02.jpg to see the difference between trimmed and untrimmed data.
Of course, if you have a solution for this, feel free to submit your version of this code which can circumvent this issue! I'm definitely interested in improving my code :)

File name format: material_load(in units of mN)_measurement#_data.txt. Example: fsilica_100mN_01_data.txt. 
With this format, you can analyze several different loads at once for a given material. If you're working with several materials, you can also store that in the software and refer to it pretty easily. But all of this hinges on this assumed file format.


