###extractMCNP.py
###Extract MCNP data from CSV file and save as TH2D in ROOT file

##File defining info
import ROOT
column = 4 #columns: 2= p>564MeV; 4=p<564MeV; 10=total
if column == 10: filename="MCNP_Pencil_570MeV_Total"
elif column == 2: filename="MCNP_Pencil_570MeV_Primary"
elif column == 4: filename="MCNP_Pencil_570MeV_Secondary"

###Load in MCNPX
##Function to read in particle densities as a 2D map:
def decode2DMap(path,file,lenx,leny,column):
    import csv
    from numpy import zeros
    Img = zeros((leny,lenx))
    Error = zeros((leny,lenx))
    j=0
    i=0
    if file[-1] != "v":
        print(file)
        file+=".csv"
        print(path+file)
    with open(path+file,mode='r') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        next(csv_reader)
        next(csv_reader)
        for row in csv_reader:
            #print(j,i,column,row)
            Img[j,i] = float(row[column])
            Error[j,i] = float(row[column+1])
            #print(j,i,Img[j,i])
            i+=1
            if i == lenx:
                i=0
                j+=1
            if j==leny+1:
                print(j,i,Img[j,i])
                print("ending",j)
                break
    csv_file.close()
    #print(Img.shape)
    return Img,Error
#Settings to read in
path = "/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/MCNPX/"
file = "MCNP_ProtonDensity_570MeV_10June"
picpath= "/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/Pictures/"
if column == 10:
    energyInfo = " Full Energy" #= column 10
    engFile = "Total"
elif column == 2:
    energyInfo = " E > 564MeV" # = column 2
    engFile = "Primary"
elif column == 4:
    energyInfo = " E < 564MeV"
    engFile = "Secondary"
lenx  = 77 #indices
lowx  = -76  #[mm]
highx = 76   #[mm]
lowy  = -200 #[mm]
highy = 200  #[mm]
leny  = 201 #indices
##Call decoder to 2D Map function
Img,Error = decode2DMap(path,file,lenx,leny,column)

###Convert MCNPX 2D map to ROOT
##Make TH2D histogram to read the 2D Map into, same size as 2D Map
MCNP_BERT = ROOT.TH2D("MCBERT_TargetProtons","Scattered Pencil Beam at Target;Horizontal [mm];Vertical [mm]",lenx,lowx,highx,leny,lowy,highy)
#Help from Kyrre's converter function and https://root-forum.cern.ch/t/create-and-fill-a-th2f-using-arrays-for-binning/27161/2
##Get each bin and set bin content and error
for xi in range(lenx):
    for yi in range(leny):
        #bx = hIn.GetBin(xi+1,yi+1)
        #MCNP_BERT.GetBin(yi,xi) = Img[yi,xi]
        bin = MCNP_BERT.GetBin(xi+1,yi+1)
        MCNP_BERT.SetBinContent(bin,Img[yi,xi])
        MCNP_BERT.SetBinError(bin,Error[yi,xi])
        #if xi == 38 and yi == 100:
            #print(Img[yi,xi],MCNP_BERT.GetBinContent(MCNP_BERT.GetBin(xi+1,yi+1)))

myFile = ROOT.TFile.Open("/scratch2/ericdf/PBWScatter/"+filename+".root","RECREATE")
myFile.WriteObject(MCNP_BERT,"MCBERT_TargetProtons")

c1 = ROOT.TCanvas("MCNP Test","MCNP Test",0,0,100*8,100*8)
newFile = ROOT.TFile("/scratch2/ericdf/PBWScatter/"+filename+".root")
test = newFile.Get("MCBERT_TargetProtons").Clone("test")
test.Draw()
c1.Draw()
c1.Print("/scratch2/ericdf/PBWScatter/"+filename+".png")
