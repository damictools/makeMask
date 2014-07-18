#include <string.h>
#include <stdio.h>
#include "fitsio.h"

#include <iostream>
#include <sstream>

#include <sys/time.h>
#include <time.h>
#include <inttypes.h>
#include <fstream>
#include <unistd.h>
#include <getopt.h>    /* for getopt_long; standard getopt is in unistd.h */
#include <vector>
#include <algorithm>
#include <ctime>
#include <climits>
#include <cmath>
#include <iomanip>

#include "globalConstants.h"

using namespace std;

int deleteFile(const char *fileName){
  cout << yellow;
  cout << "Will overwrite: " << fileName << endl << endl;
  cout << normal;
  return unlink(fileName);
}

bool fileExist(const char *fileName){
  ifstream in(fileName,ios::in);
  
  if(in.fail()){
    //cout <<"\nError reading file: " << fileName <<"\nThe file doesn't exist!\n\n";
    in.close();
    return false;
  }
  
  in.close();
  return true;
}

/*========================================================
  ASCII progress bar
==========================================================*/
void showProgress(unsigned int currEvent, unsigned int nEvent) {

  const int nProgWidth=50;

  if ( currEvent != 0 ) {
    for ( int i=0;i<nProgWidth+8;i++)
      cout << "\b";
  }

  double percent = (double) currEvent/ (double) nEvent;
  int nBars = (int) ( percent*nProgWidth );

  cout << " |";
  for ( int i=0;i<nBars-1;i++)
    cout << "=";
  if ( nBars>0 )
    cout << ">";
  for ( int i=nBars;i<nProgWidth;i++)
    cout << " ";
  cout << "| " << setw(3) << (int) (percent*100.) << "%";
  cout << flush;

}


void printCopyHelp(const char *exeName, bool printFullHelp=false){
  
  if(printFullHelp){
    cout << bold;
    cout << endl;
    cout << "This program computes the mask image from a MAD input fits file.\n";
    cout << "It handles all the available HDUs. The HDUs in the output fit file\n";
    cout << "will be 8bits and the pixels will be 0 for good ones and 8 for bad ones.\n";
    cout << "Good pixel deffinition: |pix-MEDIAN| < madCut*consistencyConst*MAD.\n";
    cout << "With consistencyConst = " << kConsistencyConstant << " so the madCut is in SD units.\n";
    cout << "\n";
    cout << normal;
  }
  cout << "==========================================================================\n";
  cout << yellow;
  cout << "\nUsage:\n";
  cout << "  "   << exeName << " <input file> -o <output filename> -c <madCut>\n";
  cout << "\nOptions:\n";
  cout << "  -q for quiet (no screen output)\n";
  cout << "  -s <HDU number> for processing a single HDU\n\n";
  cout << blue;
  cout << "For any problems or bugs contact Javier Tiffenberg <javiert@fnal.gov>\n\n";
  cout << normal;
  cout << "==========================================================================\n\n";
}

string bitpix2TypeName(int bitpix){
  
  string typeName;
  switch(bitpix) {
      case BYTE_IMG:
          typeName = "BYTE(8 bits)";
          break;
      case SHORT_IMG:
          typeName = "SHORT(16 bits)";
          break;
      case LONG_IMG:
          typeName = "INT(32 bits)";
          break;
      case FLOAT_IMG:
          typeName = "FLOAT(32 bits)";
          break;
      case DOUBLE_IMG:
          typeName = "DOUBLE(64 bits)";
          break;
      default:
          typeName = "UNKNOWN";
  }
  return typeName;
}



int copyStructure(const string &inFile, const char *outF, const int singleHdu){
    
  fitsfile  *outfptr; /* FITS file pointers defined in fitsio.h */
  fitsfile *infptr;   /* FITS file pointers defined in fitsio.h */
  
  int status = 0;
  int single = 0;
  
  int hdutype, bitpix, naxis = 0, nkeys;
  int nhdu = 0;
  long naxes[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1};
  long totpix = 0;
  //char card[81];
  
  ostringstream fileStructSummary;
  
  const char* inF = inFile.c_str();
  fits_open_file(&infptr, inF, READONLY, &status); /* Open the input file */
  if (status != 0) return(status);
  
  fits_get_num_hdus(infptr, &nhdu, &status);
  
  if (singleHdu>0){
    if(singleHdu > nhdu){
      fits_close_file(infptr,  &status);
      cerr << red << "\nError: the file does not have the required HDU!\n\n" << normal;
      return -1000;
    }
    single = 1; /* Copy only a single HDU if a specific extension was given */
  }
  
  if(single)
    fileStructSummary << "The output file will contain HDU" << singleHdu << " of " << nhdu << " availables in the input files."<< endl;
  else
    fileStructSummary << "The output file will contain all the " << nhdu << " HDUs availables in the input files.\n";
    
  
  fits_create_file(&outfptr, outF, &status);/* Create the output file */
  if (status != 0) return(status);
  
  
  fileStructSummary << bold << "HDU   hdutype  #Axis  #Cols  #Rows   IN_datatype      OUT_datatype\n" << normal;
// HDU  hdutype #Axis #Cols #Rows datatype  
  for (int n=1; n<=nhdu; ++n)  /* Main loop through each extension */
  {

    if (single) n = singleHdu;
    
    /* get image dimensions and total number of pixels in image */
    fits_movabs_hdu(infptr, n, &hdutype, &status);
    for (int i = 0; i < 9; ++i) naxes[i] = 1;
    fits_get_img_param(infptr, 9, &bitpix, &naxis, naxes, &status);
    totpix = naxes[0] * naxes[1] * naxes[2] * naxes[3] * naxes[4] * naxes[5] * naxes[6] * naxes[7] * naxes[8];
    fileStructSummary  << setw (3) << n << "  "  << setw (8) << hdutype << "  " << setw (5) << naxis << "  " << setw (5) << naxes[0] << "  " << setw (5) << naxes[1] << "  " << setw (15) << bitpix2TypeName(bitpix) << "  " << setw (15) << bitpix2TypeName( 8 );
    
    //continue;
    if (hdutype != IMAGE_HDU || naxis == 0 || totpix == 0){
      /* just copy tables and null images */
      fits_copy_hdu(infptr, outfptr, 0, &status);
      if (status != 0) return(status);
      fileStructSummary << magenta << "<- Not an image HDU" << normal;
    }
    else{
      
      fits_create_img(outfptr, BYTE_IMG, naxis, naxes, &status);
      if (status != 0) return(status);

      /* copy the relevant keywords */
      fits_get_hdrspace(infptr, &nkeys, NULL, &status); 
      for (int i = 1; i <= nkeys; ++i) {
	char card[FLEN_CARD];
        fits_read_record(infptr, i, card, &status);
	if(strncmp (card, "TRIMSEC", 7) == 0) continue;
	if(strncmp (card, "DATASEC", 7) == 0) continue;
	if(strncmp (card, "BZERO", 5) == 0) continue;
        if (fits_get_keyclass(card) > TYP_CMPRS_KEY) fits_write_record(outfptr, card, &status);
      }
      
    }
    fileStructSummary << endl;
    
    /* quit if only copying a single HDU */
    if (single) break;
  }
  fits_close_file(infptr, &status);
  fits_close_file(outfptr,  &status);

  if(gVerbosity){
    cout << bold << "Files structure summary:\n" << normal;
    cout << fileStructSummary.str();
    cout << green << "Structure copied.\n\n" << normal;
  }
  return status;
}

/* Compute the median image and write it to an existing output file 
 * that has the right structure (created by copyStructure function */
int computeImage(const string inFile, const char *outF, const int singleHdu, const float madCut){
  
  
  const float sigmaEq = kConsistencyConstant*madCut;
  
  int status = 0;
  double nulval = 0.;
  int anynul = 0;
  int single = 0;
  
  int nhdu = 0;
  
  /* Open the output file */
  fitsfile  *outfptr; /* FITS file pointers defined in fitsio.h */
  fits_open_file(&outfptr, outF, READWRITE, &status);
  if (status != 0) return(status);
  
  fits_get_num_hdus(outfptr, &nhdu, &status);
  
  if (singleHdu>0){
    single = 1; /* Copy only a single HDU if a specific extension was given */
  }
  
  fitsfile  *infptr; /* FITS file pointers defined in fitsio.h */
  const char* inF = inFile.c_str();
  fits_open_file(&infptr, inF, READONLY, &status); /* Open the input file */
  if (status != 0) return(status);
  
  
  ostringstream maskSummary;
  maskSummary << bold << "HDU   Median    MAD  #BadPix  %BadPix\n" << normal;
  
  for (int n=1; n<=nhdu; ++n)  /* Main loop through each extension */
  {
    int nOut = n;
    if (single){
      n = singleHdu;
      nOut = 1;
    }
    const int nHDUsToProcess = (single>0)? 1 : nhdu;
//     if(gVerbosity){
//       cout << bold << "\rProcessing HDU: " << n << normal << flush;
//       showProgress(0,1);
//     }
    
    /* get output image dimensions and total number of pixels in image */
    int hdutype, bitpix, bytepix, naxis = 0;
    long naxes[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1};
    fits_movabs_hdu(outfptr, nOut, &hdutype, &status);
    for (int i = 0; i < 9; ++i) naxes[i] = 1;
    fits_get_img_param(outfptr, 9, &bitpix, &naxis, naxes, &status);
    long totpix = naxes[0] * naxes[1];
    
    /* Don't try to process data if the hdu is empty */    
    if (hdutype != IMAGE_HDU || naxis == 0 || totpix == 0){
      if(single) break;
      continue;
    }
    
    bytepix = abs(bitpix) / 8;
    if(bytepix!=1) return -1000;
    
    char* outArray = new char[totpix];
    
    for(int i=0;i<totpix;++i) outArray[i] = 0;
    
    
    fits_movabs_hdu(infptr, n, &hdutype, &status);
    
    
    if(gVerbosity){
      if(single) showProgress(0,3*nHDUsToProcess);
      else showProgress((n-1)*3+0,3*nHDUsToProcess);
    }
      
    /* Open the input file */

    
    long fpixel[2]={1,1};
    long lpixel[2]={naxes[0],naxes[1]};
    
    int nLines = (lpixel[1]-fpixel[1]+1);    
    const int nCol = (lpixel[0]-fpixel[0]+1); 
    
    
    long inc[2]={1,1};
    const long npix =  nCol*nLines;
    double* sArray = new double[npix];
    /* Read the images as doubles, regardless of actual datatype. */
    fits_read_subset(infptr, TDOUBLE, fpixel, lpixel, inc, &nulval, sArray, &anynul, &status);
    
    double* sAuxArray = new double[npix];
    std::copy(sArray,sArray+npix, sAuxArray);
    sort(sAuxArray,sAuxArray+npix);
    const double sMedian = sAuxArray[npix/2];
    
    if(gVerbosity){
      if(single) showProgress(1,3*nHDUsToProcess);
      else showProgress((n-1)*3+1,3*nHDUsToProcess);
    }
    
    for(int p=0;p<npix;++p){
      sAuxArray[p] = fabs(sArray[p] - sMedian);
    }
    sort(sAuxArray,sAuxArray+npix);
    const double sMad = sAuxArray[npix/2];
    
   
    long nMaskedPix=0;
    for(int p=0;p<npix;++p){
      if( fabs(sArray[p]-sMedian) > sigmaEq*sMad || sArray[p]==0){ //Mask pixel
          outArray[p] =  128;
	  ++nMaskedPix;
	}
    }
    
    maskSummary << setw (3) << n << " " << fixed << setprecision(2) << setw (8) << sMedian << " " << setw (6) << sMad << " " << setw (8) << nMaskedPix << " " << setw (8) << nMaskedPix*100.0/npix << endl;
    
    delete[] sArray;
    
    
    
    if(gVerbosity){
      if(single) showProgress(2,3*nHDUsToProcess);
      else showProgress((n-1)*3+2,3*nHDUsToProcess);
    }

    fits_write_img(outfptr, TBYTE, 1, totpix, outArray, &status);
    

    delete[] outArray;
    

    if(gVerbosity){
      if(single) showProgress(3,3*nHDUsToProcess);
      else showProgress((n-1)*3+3,3*nHDUsToProcess);
    }
    
    
    /* quit if only copying a single HDU */
    if (single) break;
  }
  
  /* Close the output file */
  fits_close_file(outfptr,  &status);
  fits_close_file(infptr,  &status);
  if(gVerbosity){
    showProgress(1,1);
  }
  
  if(gVerbosity){
    cout << bold << "\n\nMask summary:\n" << normal;
    cout << maskSummary.str();
  }
  return status;
}


void checkArch(){
  if(sizeof(float)*CHAR_BIT!=32 || sizeof(double)*CHAR_BIT!=64){
    cout << red;
    cout << "\n ========================================================================================\n";
    cout << "   WARNING: the size of the float and double variables is non-standard in this computer.\n";
    cout << "   The program may malfunction or produce incorrect results\n";
    cout << " ========================================================================================\n";
    cout << normal;
  }
}

int processCommandLineArgs(const int argc, char *argv[], int &singleHdu, float &madCut, string &inFile, string &outFile){
  
  if(argc == 1) return 1;
  
  bool outFileFlag = false;
  string inListFile = "";
  int opt=0;
  while ( (opt = getopt(argc, argv, "c:o:s:qQhH?")) != -1) {
    switch (opt) {
    case 'o':
      if(!outFileFlag){
        outFile = optarg;
        outFileFlag = true;
      }
      else{
        cerr << red << "\nError, can not set more than one output file!\n\n" << normal;
        return 2;
      }
      break;
    case 's':
      if(singleHdu<0){
        singleHdu = atoi(optarg);
      }
      else{
        cerr << red << "\nError, can not set more than one HDU!\n\n"  << normal;
        return 2;
      }
      break;
    case 'c':
      if(madCut<0){
        madCut = atof(optarg);
      }
      else{
        cerr << red << "\nError, can not set more than madCut!\n\n"  << normal;
        return 2;
      }
      break;
    case 'Q':
    case 'q':
      gVerbosity = 0;
      break;
    case 'h':
    case 'H':
    default: /* '?' */
      return 1;
    }
  }
  
  if(!outFileFlag){
    cerr << red << "\nError: output filename missing.\n" << normal;
    return 2;
  }
  
  if(madCut<0){
    cerr << red << "Error: missing or invalid madCut!\n\n" << normal;
    return 1;
  }

  inFile="";
  
  if(argc-optind==0){
    cerr << red << "Error: no input file(s) provided!\n\n" << normal;
    return 1;
  }
  else if(argc-optind>1){
    cerr << red << "Error: more than one input file provided!\n\n" << normal;
    return 1;
  }
  
  inFile=argv[optind];
  if(!fileExist(inFile.c_str())){
    cout << red << "\nError reading input file: " << inFile <<"\nThe file doesn't exist!\n\n" << normal;
    return 1;
  }
  
  return 0;
}

int main(int argc, char *argv[])
{
  
  checkArch(); //Check the size of the double and float variables.
  
  time_t start,end;
  double dif;
  time (&start);
  
  string outFile;
  string inFile;
  int singleHdu=-1;
  float madCut=-1;
  
  int returnCode = processCommandLineArgs( argc, argv, singleHdu, madCut, inFile, outFile);
  if(returnCode!=0){
    if(returnCode == 1) printCopyHelp(argv[0],true);
    if(returnCode == 2) printCopyHelp(argv[0]);
    return returnCode;
  }
  
  if(gVerbosity){
    cout << bold << "\nWill read the following file:\n" << normal;
    cout << "\t" << inFile << endl;
    cout << bold << "\nThe output will be saved in the file:\n\t" << normal << outFile << endl;
    cout << bold << "\nWill use: " << cyan << kConsistencyConstant <<"x"<<  madCut <<"xMAD "<< normal << bold << "as pixel cut value." << endl;
  }
  
  /* Overwrite the output file if it already exist */
  if(fileExist(outFile.c_str())){
    cout << yellow << "\nThe output file exist. " << normal;
    deleteFile(outFile.c_str());
  }
  
  
  /* Do the actual processing */
  int status = copyStructure( inFile,  outFile.c_str(), singleHdu);
  if (status != 0){ 
    fits_report_error(stderr, status);
    return status;
  }
  
  
  status = computeImage( inFile,  outFile.c_str(), singleHdu, madCut);
  if (status != 0){ 
    fits_report_error(stderr, status);
    return status;
  }
  
  /* Report */
  time (&end);
  dif = difftime (end,start);
  if(gVerbosity) cout << green << "All done!\n" << bold << "-> It took me " << dif << " seconds to do it!\n\n" << normal;

  return status;
}
