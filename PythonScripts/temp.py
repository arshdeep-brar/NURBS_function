import csv

Data_pts = []

with open('NACA-2412.csv') as csvfile:
     readCSV = csv.reader(csvfile, delimiter=',')
     for row in readCSV:
         Data_pts.append([float(row[0]),float(row[1])])

        

