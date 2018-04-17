#if !defined(linux) && !defined(__linux) && !defined(__linux__)
#   error Sorry! Linux only code!
#endif // #if !defined(linux) && !defined(__linux) && !defined(__linux__)
#include <stdio.h>
#include <math.h>
#include <cstdlib>
#include <unistd.h>
#include <iostream>
#include <errno.h>
#include <getopt.h>
#include <apps/Common/exampleHelper.h>
#include <mvIMPACT_CPP/mvIMPACT_acquire.h>
#include <mvIMPACT_CPP/mvIMPACT_acquire_GenICam.h>

#include <fitsio.h>
#include <fitsio2.h>

void fits_report_error (int status);

using namespace std;
using namespace mvIMPACT::acquire;

#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))

/******************
globals for use in coordinating results from 4 cameras
******************/
int g_verbose=0;
int g_save_raw=0;
int g_save_rebin=0;
int g_save_mosaic=0;
int g_pLock=0;
unsigned int *g_pMosaic=NULL;
int g_decimate=0;
int g_mosaic_rebin=4;
int g_mosaicNX=0;
int g_mosaicNY=0;
int g_mosaic_tile[]={2,2};
char **g_serials;
int g_ncameras=0;
long g_expotime=200000;
char g_lvl=00;
/******************
end of globals
******************/

/******************
fits related structures
******************/

typedef struct proj_x_s {
  int   nbin;
  float *proj_x_coord;
  int   *proj_x_npix;
  int   *proj_x_hist;
  float *proj_x_density;
} proj_x;
  
typedef struct map_info_s {
  long   naxes[2];
  char   *map;
  int    n_px;
  proj_x *px;
} map_info;

typedef struct fits_s {
  char *filename;
  fitsfile *ff;
  map_info m_info;
  struct fits_s *next_fits;
} fits;

void append_fits_request (fits **f,char *name,Request *thisreq);
void destroy_fits_request (fits **f);

/******************
end of fits
******************/

/******************
roi related structures
******************/

typedef struct roi_s {
  int   camera_index;
  float rot;
  float lim[2][2]; // first index is the corner, second is the axis.
  struct roi_s *next_roi;
  // other items that are derived from the image in question and this extraction spec.
} roi;

void append_roi  (roi **r,char *roispec);
void report_rois (roi *r);

static roi *g_roih=NULL;


/******************
end of roi
******************/


/******************
routines relevant to both roi and fits
******************/

void compute_roi(roi *rp,fits *fp,int roi_ix);
void compute_roi_fits_request (roi *r,fits *f,int camera_index);

/******************
end of roi / fits
******************/

void stash_results ( Request* thisreq, Device* pDev);

#define PRESS_A_KEY				\
  getchar();

//-----------------------------------------------------------------------------
class ThreadParameter
//-----------------------------------------------------------------------------
{
  Device*         m_pDev;
  volatile bool   m_boTerminateThread;
public:
  ThreadParameter( Device* pDev ) : m_pDev( pDev ), m_boTerminateThread( false ) {}
  Device* device( void ) const
  {
    return m_pDev;
  }
  bool    terminated( void ) const
  {
    return m_boTerminateThread;
  }
  void    terminateThread( void )
  {
    m_boTerminateThread = true;
  }
};

//-----------------------------------------------------------------------------
static unsigned int thread_func( void* pData )
//-----------------------------------------------------------------------------
{
  ThreadParameter* pThreadParameter = reinterpret_cast<ThreadParameter*>( pData );
  Device *pDev = pThreadParameter->device();

  unsigned int cnt = 0;
  
  try
    {
      // this is the place to set interface etc.
      // need to use GenICam layout for BlueFox3 devices, device specific layout is deprecated!
      conditionalSetProperty(pDev->interfaceLayout,dilGenICam);
      pDev->open();
    }
  catch( const ImpactAcquireException& e )
    {
      // this e.g. might happen if the same device is already opened in another process...
      cout << "An error occurred while opening the device " << pDev->serial.read() << "(error code: " << e.getErrorCode() << "(" << e.getErrorCodeAsString() << ")). Terminating thread." << endl << "Press [ENTER] to end the application..." << endl;
      PRESS_A_KEY
        return 0;
    }
  

  GenICam::AcquisitionControl ac( pDev );
  
  cout << "for " << pDev->serial.read()  << " exposureMode (current): " << ac.exposureMode.readS() << endl;
  ac.exposureMode.writeS("Timed");
  cout << "exposureMode (new): " << ac.exposureMode.readS() << endl;

  {
    GenICam::AnalogControl agc( pDev );
    if (g_verbose) 
      cout << "for " << pDev->serial.read()  << "gainAuto (current): " << agc.gainAuto.readS() << endl;
    agc.gainAuto.writeS( "Off" );
    if (g_verbose) 
      cout << "gainAuto (new): " << agc.gainAuto.readS() << endl;
    
    if (g_verbose) 
      cout << "for " << pDev->serial.read()  << "gain (current): " << agc.gain.read() << endl;
    agc.gain.write( 1 );
    if (g_verbose) 
      cout << "gain (new): " << agc.gain.read() << endl;
  }

  if ((g_expotime>0) || (g_lvl>0)) { // for timed exposures
    cout << "for " << pDev->serial.read()  << "exposureAuto (current): " << ac.exposureAuto.readS() << endl;
    ac.exposureAuto.writeS( "Off" );
    cout << "exposureAuto (new): " << ac.exposureAuto.readS() << endl;
    
    // the following works!!
    cout << "for " << pDev->serial.read()  << "exposureTime (current): " << ac.exposureTime.read() << endl;
    ac.exposureTime.write( g_expotime );
    cout << "exposureTime (new): " << ac.exposureTime.read() << endl;
  } else { // for auto exposure
    cout << "for " << pDev->serial.read()  << "exposureAuto (current): " << ac.exposureAuto.readS() << endl;
    ac.exposureAuto.writeS("Continuous");
    cout << "exposureAuto (new): " << ac.exposureAuto.readS() << endl;
  }
  
  GenICam::ImageFormatControl ifc( pDev );
  if (g_verbose) 
    cout << "for " << pDev->serial.read()  << "pixelFormat (current): " << ifc.pixelFormat.readS() << endl;
  ifc.pixelFormat.writeS("Mono8");
  if (g_verbose) 
    cout << "pixelFormat (new): " << ifc.pixelFormat.readS() << endl;
  


  // establish access to the statistic properties
  Statistics statistics( pDev );
  // create an interface to the device found
  FunctionInterface fi( pDev );

  
  // pre-fill the capture queue. There can be more than 1 queue for some devices, but for this sample
  // we will work with the default capture queue. If a device supports more than one capture or result
  // queue, this will be stated in the manual. If nothing is mentioned about it, the device supports one
  // queue only. Request as many images as possible. If there are no more free requests 'DEV_NO_FREE_REQUEST_AVAILABLE'
  // will be returned by the driver.
  int result = DMR_NO_ERROR;
  SystemSettings ss( pDev );
  const int REQUEST_COUNT = ss.requestCount.read();
  for( int i = 0; i < REQUEST_COUNT; i++ )
    {
      result = fi.imageRequestSingle();
      if( result != DMR_NO_ERROR )
        {
	  cout << "Error while filling the request queue: " << ImpactAcquireException::getErrorCodeAsString( result ) << endl;
        }
    }

  // run thread loop
  const Request* pRequest = 0;
  const unsigned int timeout_ms = 8000;   // USB 1.1 on an embedded system needs a large timeout for the first image
  int requestNr = INVALID_ID;
  // This next comment is valid once we have a display:
  // we always have to keep at least 2 images as the display module might want to repaint the image, thus we
  // can't free it unless we have a assigned the display to a new buffer.
  int lastRequestNr = INVALID_ID;
  while( !pThreadParameter->terminated() )
    {
      // wait for results from the default capture queue
      requestNr = fi.imageRequestWaitFor( timeout_ms );
      if( fi.isRequestNrValid( requestNr ) )
        {
	  pRequest = fi.getRequest( requestNr );
	  if( pRequest->isOK() )
            {
	      ++cnt;
	      // here we can display some statistical information every 100th image
	      if( cnt % 100 == 0 )
                {
		  cout << "Info from " << pDev->serial.read()
		       << ": " << statistics.framesPerSecond.name() << ": " << statistics.framesPerSecond.readS()
		       << ", " << statistics.errorCount.name() << ": " << statistics.errorCount.readS()
		       << ", " << statistics.captureTime_s.name() << ": " << statistics.captureTime_s.readS() << endl;
                }
            }
	  else
            {
	      cout << "Error: " << pRequest->requestResult.readS() << endl;
            }

	  if( fi.isRequestNrValid( lastRequestNr ) )
	    {
	      // save as fits file
	      stash_results( fi.getRequest(lastRequestNr), pDev );
	      // the buffer is no longer needed...
	      fi.imageRequestUnlock( lastRequestNr );
	    }
	  lastRequestNr = requestNr;
	  // send a new image request into the capture queue
	  fi.imageRequestSingle();
        }
      else
        {
	  // If the error code is -2119(DEV_WAIT_FOR_REQUEST_FAILED), the documentation will provide
	  // additional information under TDMR_ERROR in the interface reference (
	  cout << "imageRequestWaitFor failed (" << requestNr << ", " << ImpactAcquireException::getErrorCodeAsString( requestNr ) << ", device " << pDev->serial.read() << ")"
	       << ", timeout value too small?" << endl;
        }
    }

  // free the last potentially locked request
  if( fi.isRequestNrValid( requestNr ) )
    {
      fi.imageRequestUnlock( requestNr );
    }
  // clear all queues
  fi.imageRequestReset( 0, 0 );
  return 0;
}

//-----------------------------------------------------------------------------
static void* liveThread( void* pData )
//-----------------------------------------------------------------------------
{
  thread_func( pData );
  return NULL;
}

//-----------------------------------------------------------------------------
int main( int argc, char* argv[] )
//-----------------------------------------------------------------------------
{
  // process commandline arguments
  int c;

  while (1) {

    static struct option long_options[]=
      {{"verbose", no_argument,       &g_verbose,      1},
       {"decimate",no_argument,       &g_decimate,     1},
       {"bin",     required_argument, 0,             'b'},
       {"expoauto",required_argument, 0,             'a'},
       {"expotime",required_argument, 0,             'x'},
       {"save_raw",no_argument,       0,             's'},
       {"mosaic",  required_argument, 0,             'm'},
       {"roi",     required_argument, 0,             'r'},
       {0,         0,                 0,             0}};

    int option_index=0;
    c=getopt_long(argc,argv,"bdsm:r:x:",long_options,&option_index);

    if (c==-1) break;

    switch (c) {
    case 0:
      if (long_options[option_index].flag != 0) break;
      printf ("option %s",long_options[option_index].name);
      if (optarg) printf(" with arg %s",optarg);
      printf ("\n");
      break;

    case 'a':
      printf("option -a (expoauto) with value `%s'\n",optarg);
      sscanf(optarg,"%c",&g_lvl);
      break;

    case 'x':
      printf("option -x (expotime) with value `%s'\n",optarg);
      sscanf(optarg,"%ld",&g_expotime);
      break;

    case 'b':
      printf("option -b (bin) with value `%s'\n",optarg);
      sscanf(optarg,"%d",&g_save_rebin);
      break;

    case 's':
      printf("option -s (save raw)\n");
      g_save_raw=1;
      break;

    case 'd':
      printf("option -d (decimate)\n");
      g_decimate=1;
      break;

    case 'm':
      printf("option -m (mosaic) with rebin value `%s'\n",optarg);
      g_save_mosaic=1;
      sscanf(optarg,"%d",&g_mosaic_rebin);
      break;

    case 'r':
      printf("option -r (ROI specification) with value `%s'\n",optarg);
      append_roi(&g_roih,optarg);
      break;

    case '?':
      // do nothing, already printed error
      break;
    default:
      abort();
      
    }
  }

  if (g_verbose)
    puts ("verbose flag is set");

  report_rois(g_roih);

  if (optind < argc) {
    printf ("non-option ARGV-elements: ");
    //      printf ("non-option ARGV-elements (these are all fits files: ");
    //      while (optind < argc) {
    //	append_fits(&fitsh,argv[optind++]);
    //      }
    putchar('\n');
  }

  // get on to device manager stuff(s)
  DeviceManager devMgr;
  const unsigned int devCnt = devMgr.deviceCount();
  if( devCnt == 0 )
    {
      cout << "No MATRIX VISION device found! Unable to continue!" << endl;
      return 0;
    }

  pthread_t* pHandles = new pthread_t[devCnt];
  pthread_attr_t* pAttrs = new pthread_attr_t[devCnt];
  // store all device infos in a vector
  // and start the execution of a 'live' thread for each device.
  vector<ThreadParameter*> threadParams;
  for( unsigned int i = 0; i < devCnt; i++ )
    {
      threadParams.push_back( new ThreadParameter( devMgr[i] ) );
      cout << devMgr[i]->family.read() << "(" << devMgr[i]->serial.read() << ")" << endl;
    }
  
  { // set up some globals for coordinated output
    g_ncameras=devCnt;
    g_serials=(char**)malloc(g_ncameras*sizeof(char*));
    for (int j=0;j<g_ncameras;j++) {
      g_serials[j]=(char*)malloc(1024*sizeof(char));
    }
    g_mosaic_tile[1] = ceil(g_ncameras / g_mosaic_tile[0]);
    for( unsigned int i = 0; i < devCnt; i++ )
      {
	// populate globals for use later
	sprintf(g_serials[i],devMgr[i]->serial.read().c_str());
      }
  }

  
  // start live threads
  for( unsigned int j = 0; j < devCnt; j++ )
    {
      pthread_attr_init( &pAttrs[j] );
      // you can set the stack size like this: pthread_attr_setstacksize (&pAttrs[j], 1024*1024);
      pthread_create( &pHandles[j], &pAttrs[j], liveThread, ( void* )threadParams[j] );
    }

  // now all threads will start running...
  cout << "Press return to end the acquisition( the initialisation of the devices might take some time )" << endl;
  PRESS_A_KEY

    // stop all threads again
    cout << "Terminating live threads..." << endl;
  size_t vSize = threadParams.size();
  for( unsigned int k = 0; k < vSize; k++ )
    {
      cout << "Terminating thread " << k << "." << endl;
      threadParams[k]->terminateThread();
    }

  // wait until each live thread has terminated.
  for( unsigned int j = 0; j < devCnt; j++ )
    {
      cout << "Waiting for thread " << j << " to terminate." << endl;
      pthread_join ( pHandles[j], NULL );
    }
  cout << "All capture threads terminated." << endl;

  // free resources
  for( unsigned int l = 0; l < vSize; l++ )
    {
      delete threadParams[l];
      cout << "thread parameter " << l << " removed." << endl;
    }
  delete [] pHandles;
  cout << "All thread handles removed." << endl;

  delete [] pAttrs;
  cout << "All thread atributes removed." << endl;
  return 0;
}

void stash_results (Request* thisreq, Device *pDev) {

  if (thisreq->isOK()) {
    unsigned int ts = thisreq->infoTimeStamp_us.read();
    unsigned int xt = thisreq->infoExposeTime_us.read();
    string ser = pDev->serial.read();
    //    unsigned int frame_ix = thisreq->infoFrameID.read();
    const char* pData = (char*) thisreq->imageData.read();
    unsigned int iWidth = thisreq->imageWidthTotal.read();
    unsigned int iHeight = thisreq->imageHeightTotal.read();
    unsigned int npix=iWidth*iHeight;

    if (g_save_mosaic) {
      // figure out which sensor this is (by serial) to determine its place
      int camera_index=g_ncameras;
      while (camera_index-- && strcmp(g_serials[camera_index],pDev->serial.read().c_str()));
      if (g_verbose)
	cout << "got camera index: " << camera_index << endl;
      // initialize if necessary, using g_pLock etc.
      while (g_pLock) {
	// do nothing until lock is released
	if (g_verbose)
	  cout << " waiting .. (camera_index=" << camera_index << ")";
	usleep(30);
      }
      if (g_pMosaic==NULL) {
	g_pLock=1; // lock to prevent multiple allocation
	// need to allocate for the mosaic image
	if (g_verbose)
	  cout << "locking..";
	g_mosaicNX=g_mosaic_tile[0]*(iWidth/g_mosaic_rebin);
	g_mosaicNY=g_mosaic_tile[1]*(iHeight/g_mosaic_rebin);
	g_pMosaic=(unsigned int*)malloc(g_mosaicNX*g_mosaicNY*sizeof(unsigned int));
	g_pLock=0;
	if (g_verbose)
	  cout << "unlocked!" << endl;
      }
      // proceed to rebin this image and place into g_pMosaic
      {
	int startx=(camera_index%g_mosaic_tile[0])*(iWidth/g_mosaic_rebin);
	int stopx=startx+(iWidth/g_mosaic_rebin);
	int starty=(camera_index/g_mosaic_tile[0])*(iHeight/g_mosaic_rebin);
	int stopy=starty+(iHeight/g_mosaic_rebin);

	if (g_verbose)
	  cout << "for camera index " << camera_index << 
	    " have (startx,stopx,starty,stopy) = (" << 
	    startx << "," << stopx << "," << starty << "," << stopy << ")" << endl;
	
	while (g_pLock) {
	  // do nothing until lock is released
	  if (g_verbose)
	    cout << " waiting .. (camera_index=" << camera_index << ")";
	  usleep(30);
	}
	g_pLock=1;
	// clear and then populate with current image
	for (int y=starty;y<stopy;y++) {
	  for (int x=startx;x<stopx;x++) {
	    g_pMosaic[y*g_mosaicNX+x]=0;
	  }
	}
	unsigned long i;
	i=iHeight*iWidth;
	while (i--) {
	  int x=startx+(i%iWidth)/g_mosaic_rebin;
	  int y=starty+(i/iWidth)/g_mosaic_rebin;
	  g_pMosaic[y*g_mosaicNX+x] += (0xFF & ((char*)(pData))[i]);
	}
	// save the mosaic if appropriate (do this only for a specific camera index case
	if (camera_index==0) {
	  // make a copy of the data and unlock for other camera threads to proceed
	  unsigned int* mosaic_copy=(unsigned int*)malloc(g_mosaicNX*g_mosaicNY*sizeof(unsigned int));
	  memcpy(mosaic_copy,g_pMosaic,g_mosaicNX*g_mosaicNY*sizeof(unsigned int));
	  g_pLock=0;
	  // unfinished business of saving the fits file and freeing up the copy image.
	  int status=0;
	  char outfile[1024];
	  fitsfile *ff;
	  sprintf(outfile,"dat/BF_%s_%u_%u.fits","mosaic",ts,xt);
	  if (1 || g_verbose)
	    cout << " (" << outfile << ") " << endl;
	  
	  fits_create_file(&ff,outfile,&status);
	  if (status) fits_report_error(status);
	  long naxis=2;
	  long naxes[]={(long)g_mosaicNX,(long)g_mosaicNY};
	  fits_create_img(ff,USHORT_IMG,naxis,naxes,&status);
	  if (status) fits_report_error(status);
	  fits_write_img(ff,TUINT,1,naxes[0]*naxes[1],(void*)mosaic_copy,&status);
	  if (status) fits_report_error(status);
	  fits_close_file(ff,&status);
	  if (status) fits_report_error(status);
	  free(mosaic_copy);
	} else {
	  // if other camera_index, done with contributing to g_pMosaic
	  g_pLock=0;
	}
      }
    }

    if (1) {
      unsigned long i;
      unsigned long tot_pix=0L,tot_sig=0L;
      int hist[256],chist[256];
      for (i=0L;i<256L;i++) {
	hist[i]=0;
      }
      for (i=0L;i<npix;i++) {
	unsigned char dat=((char*)(pData))[i];
	hist[dat]++;
	tot_sig+=pData[i];
	tot_pix++;
      }
      //      float avg_sig = tot_sig/((float)(tot_pix));
      // get quantile levels for adjusting exposure times
      const int nlvls=3;
      int ilvl;
      float lvls[]={0.1,0.5,0.9};
      float sig_lvls[nlvls];

      for (ilvl=0;ilvl<nlvls;ilvl++) {
	sig_lvls[ilvl]=-1.0;
      }

      for (i=0L;i<256L;i++) {
	chist[i]=(i==0)?hist[i]:chist[i-1]+hist[i];
      }

      for (ilvl=0;ilvl<nlvls;ilvl++) {
	float lvl=lvls[ilvl]*npix;
	sig_lvls[ilvl]=0;
	for (i=1L;i<256L;i++) {
	  if ((chist[i]-lvl)*(lvl-chist[i-1])>0) {
	    // cumulative distribution is just crossing lvl..
	    sig_lvls[ilvl]=(i-1)+(lvl-chist[i-1])/((float)(chist[i]-chist[i-1]));
	    goto next_lvl;
	  }
	}
      next_lvl:
	continue;
      }

      if (g_lvl>0) {
	// try out a new exposure time to aim better..
	int aim_lvl[]={(int)10,(int)100,(int)240};
	int this_lvl[]={(int)sig_lvls[0],(int)sig_lvls[1],(int)sig_lvls[2]};
	float scalar;

	int which=-1;
	which=1;

	if ((which-0)*(2-which)>=0) {
	  scalar=aim_lvl[which]/(float)(this_lvl[which]+1.0e-4);
	  scalar=g_lvl/(float)(this_lvl[which]+1.0e-4);
	} else {
	  scalar=1.0;
	}

	GenICam::AcquisitionControl ac( pDev );
	float new_exposure_time = scalar * ac.exposureTime.read();
	if (g_verbose)
	  cout << "adjust exposure time (current) " << ac.exposureTime.read() << " to scale by " << scalar << endl;
	if (isinf(new_exposure_time)) {
	  if (g_verbose)
	    cout << "not scaling.. " << endl;
	} else {
	  if ((new_exposure_time-1.0e4)*(1.0e6-new_exposure_time)>0) {
	    if (g_verbose)
	      cout << "scaling.. new exposure time " << new_exposure_time << endl;
	  } else {
	    if (new_exposure_time > 1.0e6)
	      new_exposure_time=1.0e6;
	    if (new_exposure_time < 1.0e4)
	      new_exposure_time=1.0e4;
	    if (g_verbose)
	      cout << " instead using exposure time " << new_exposure_time << endl;
	  }
	  ac.exposureTime.write(new_exposure_time);
	}
      }

      if (g_verbose)
	cout << "chist total: " << chist[255] << " -- ";
      for (ilvl=0;ilvl<nlvls;ilvl++) {
	if (g_verbose)
	  cout << lvls[ilvl] << ": " << sig_lvls[ilvl] << "; ";
      }
      if (0) { // this was just for debugging
	// and save the histogram for now
	char histfile[1024];
	sprintf(histfile,"dat/BF_%s_%u.hist",ser.c_str(),ts);
	FILE *fp;
	fp=fopen(histfile,"w");
	for (i=1L;i<256L;i++) {
	  fprintf(fp,"%d %g\n",(int)(i),chist[i]/(float)npix);
	}      
	fclose(fp);
      }
    }
    
    char outfile[1024];
    sprintf(outfile,"dat/BF_%s_%u_%u.fits",ser.c_str(),ts,xt);

    if (g_save_raw) {
      int status=0;
      fitsfile *ff;
      // fills disk space quickly and can maintain ~ 33 fps * 5Mbyte on 4 cameras simultaneously without impacting frame rate
      if (g_verbose)
	cout << " (" << outfile << ") " << endl;
      fits_create_file(&ff,outfile,&status);
      if (status) fits_report_error(status);
      long naxis=2;
      long naxes[]={(long)iWidth,(long)iHeight};
      fits_create_img(ff,BYTE_IMG,naxis,naxes,&status);
      if (status) fits_report_error(status);
      fits_write_img(ff,TBYTE,1,naxes[0]*naxes[1],(void*)pData,&status);
      if (status) fits_report_error(status);
      fits_close_file(ff,&status);
      if (status) fits_report_error(status);
    }

    if (g_save_rebin!=0) { 
      // save another representation of the image, 
      // this time decimated by g_save_rebin x g_save_rebin for smaller size

      char outfile_rebin[1024];

      if (g_decimate) {
	sprintf(outfile_rebin,"dat/BF_%s_dm%02d_%u_%u.fits",
		ser.c_str(),g_save_rebin,ts,xt);
      } else {
	sprintf(outfile_rebin,"dat/BF_%s_rb%02d_%u_%u.fits",
		ser.c_str(),g_save_rebin,ts,xt);
      }

      int status=0;
      fitsfile *ff;

      if (g_verbose)
	cout << " (" << outfile_rebin << ") " << endl;

      fits_create_file(&ff,outfile_rebin,&status);
      if (status) fits_report_error(status);

      long naxis=2;
      long naxes[]={(long)floor(iWidth/g_save_rebin),(long)floor(iHeight/g_save_rebin)};

      if (g_decimate) 
	{
	  int ix=naxes[0]*naxes[1];
	  int i,j;
	  char *pData_rebin=
	    (char*)malloc(naxes[0]*naxes[1]*sizeof(char));
	  while (ix--) {
	    i=(ix%naxes[0])*g_save_rebin;
	    j=(ix/naxes[0])*g_save_rebin;
	    pData_rebin[ix]=((char*)(pData))[j*iWidth+i];
	  }
	  fits_create_img(ff,BYTE_IMG,naxis,naxes,&status);
	  if (status) fits_report_error(status);
	  fits_write_img(ff,TBYTE,1,naxes[0]*naxes[1],
			 (void*)pData_rebin,&status);
	  if (status) fits_report_error(status);
	  fits_close_file(ff,&status);
	  if (status) fits_report_error(status);
	  free(pData_rebin);
	}
      if (!g_decimate)
	{
	  unsigned short *pData_rebin=
	    (unsigned short*)malloc(naxes[0]*naxes[1]*
				    sizeof(unsigned short));
	  int ix=naxes[0]*naxes[1];
	  while (ix--) pData_rebin[ix]=0;
	  ix=iWidth*iHeight;
	  while (ix--) {
	    int i,j;
	    i=(ix%iWidth)/g_save_rebin;
	    j=(ix/iWidth)/g_save_rebin;
	    //	    cout << "orig_x,y = " << ix%iWidth << "," << ix/iWidth << "rebin_x,y = " << i << "," << j << endl;
	    pData_rebin[j*naxes[0]+i] += (0xFF & ((char*)(pData))[ix]);
	  }
	  fits_create_img(ff,USHORT_IMG,naxis,naxes,&status);
	  if (status) fits_report_error(status);
	  fits_write_img(ff,TUSHORT,1,naxes[0]*naxes[1],
			 (void*)pData_rebin,&status);
	  if (status) fits_report_error(status);
	  fits_close_file(ff,&status);
	  if (status) fits_report_error(status);
	  free(pData_rebin);
	}
    }
    // finally compute rois on pData
    {
      // append_fits_request gets its information & buffer from memory, not from the 
      // output file just saved (if in fact it was saved). name of fits file (outfile)
      // is used primarily for display purposes.
      // "fitsh" a single element fits list (be sure to deallocate after use)
      int camera_index=g_ncameras;
      while (camera_index-- && strcmp(g_serials[camera_index],pDev->serial.read().c_str()));

      fits *fitsh=NULL;
      append_fits_request(&fitsh,outfile,thisreq); 
      compute_roi_fits_request(g_roih,fitsh,camera_index);
      destroy_fits_request(&fitsh);
    }
  }
}


void compute_roi_fits_request(roi *r,fits *f,int cix) {
  // allocate and deallocate anything within the roi & fits chains. 
  fits *fp;
  roi  *rp;
  int n_roi=0,n_fits=0;

  if (g_verbose)
    cout << "for camera index " << cix << " r= " << r << " f= " << f << endl;
  rp=r;

  while (rp!=NULL) {
    int thiscamera=rp->camera_index;
    if (cix == thiscamera) n_roi++; 
    rp=rp->next_roi;
  }

  if (n_roi==0) return;
  
  fp=f;
  while (fp!=NULL) {
    int roi_ix;
    n_fits++;
    fp->m_info.n_px=n_roi;
    if (fp->m_info.px == NULL)
      fp->m_info.px=(proj_x*)malloc(fp->m_info.n_px*sizeof(proj_x));

    rp=r;
    roi_ix=0;
    while (rp!=NULL) {
      // important - apply the same filter as above here
      int thiscamera=rp->camera_index;
      if (cix == thiscamera) {
	compute_roi(rp,fp,roi_ix);
	roi_ix++;
      }
      rp=rp->next_roi;
    }
    fp=fp->next_fits;
  }
  if (g_verbose)
    cout << "n_roi " << n_roi << " n_fits " << n_fits << endl;
  // do work
  // BLAH BLAH BLAH
  // deallocate
  fp=f;
  while (fp!=NULL) {
    if (fp->m_info.n_px==0) {
      // do nothing
    } else {
      int roi_ix;
      if (g_verbose)
	cout << "this " << cix << " has n_px = " << fp->m_info.n_px << endl;
      for (roi_ix=0;roi_ix<fp->m_info.n_px;roi_ix++) {
	// do nothing
      }
      free(fp->m_info.px);
      fp->m_info.px=NULL;
    }
    fp=fp->next_fits;
  }
}

void append_roi (roi **r,char *roispec) {
  roi *new_roi=(roi*)malloc(sizeof(roi));

  int n_parse=sscanf(roispec,"%d,%f:%f,%f:%f,%f",
		     &new_roi->camera_index,
		     &new_roi->lim[0][0],&new_roi->lim[1][0],
		     &new_roi->lim[0][1],&new_roi->lim[1][1],
		     &new_roi->rot);
  if (n_parse!=6) {
    cout << "didn't parse the expected # of fields from roispec: " <<roispec<<endl;
    cout << "expected 6, got " << n_parse << endl << "exiting.." << endl;
    exit(1);
  } else {
    if (g_verbose)
      cout << "parsed " << n_parse << " fields from roispec " << roispec << endl;
  }

  if (0) { // the following puts the current roi spec on the head of the roi list:
    new_roi->next_roi=*r;
    *r=new_roi;
  } else { // append the current roi spec, keeping the head as is.
    roi **tmprh=r;
    while (*tmprh != NULL) tmprh=&(*tmprh)->next_roi;
    *tmprh = new_roi;
  }
}

void report_rois (roi *r) {
  printf ("current contents of roi:\n");
  roi *this_roi=r;
  while (this_roi != NULL) {
    printf ("camera index = %d\n",this_roi->camera_index);
    printf ("%f < x < %f\n",this_roi->lim[0][0],this_roi->lim[1][0]);
    printf ("%f < y < %f\n",this_roi->lim[0][1],this_roi->lim[1][1]);
    printf ("theta = %f\n",this_roi->rot);
    this_roi=this_roi->next_roi;
  } 
}

void destroy_fits_request (fits **f) {
  while (*f != NULL) {
    fits *fp=*f;
    *f=(*f)->next_fits;
    free(fp);
  }
}

void append_fits_request (fits **f,char *name,Request *thisreq) {
  fits *new_fits=(fits*)malloc(sizeof(fits));
  new_fits->filename=name;
  new_fits->m_info.px=NULL;
  // copy the map over from the request (already in memory)
  new_fits->m_info.naxes[0]=thisreq->imageWidthTotal.read();
  new_fits->m_info.naxes[1]=thisreq->imageHeightTotal.read();
  new_fits->m_info.map=(char*)thisreq->imageData.read(); // isn't malloc, don't free
  new_fits->next_fits=NULL;
  if (0) { // the following puts the current roi spec on the head of the roi list:
    new_fits->next_fits=*f;
    *f=new_fits;
  } else { // append the current roi spec, keeping the head as is.
    fits **tmpfh=f;
    while (*tmpfh != NULL) tmpfh=&(*tmpfh)->next_fits;
    *tmpfh = new_fits;
  }
  // and now report the current fits list
  fits *this_fits=*f;
  if (g_verbose)
    cout << "this fits: ";
  while (this_fits != NULL) {
    if (g_verbose)
      cout << this_fits->filename;
    this_fits=this_fits->next_fits;
  }
  if (g_verbose)  
    cout << " (end)" << endl;
}

void compute_roi(roi *rp,fits *fp,int roi_ix) {
  float x0=rp->lim[0][0];    float x1=rp->lim[1][0];
  float y0=rp->lim[0][1];    float y1=rp->lim[1][1];
  float theta=rp->rot;
  // avoid boundaries that violate limits
  int x_s = (int)((x0<0)?0:x0);
  int x_f = (int)((x1<fp->m_info.naxes[0]-1)?x1:fp->m_info.naxes[0]-1);
  int y_s = (int)((y0<0)?0:y0);
  int y_f = (int)((y1<fp->m_info.naxes[1]-1)?y1:fp->m_info.naxes[1]-1);
  // center coordinates for this roi
  float x_c=(x_s+x_f)/2.0;
  float y_c=(y_s+y_f)/2.0;
  // actually use the intended center, that's where this needs to rotate about
  x_c=(x0+x1)/2.0;
  y_c=(y0+y1)/2.0;
  // find the limits (in pixels) when the roi is rotated as prescribed
  float corners[][2]={{(float)(x_s-0.5),(float)(y_s-0.5)},
		      {(float)(x_f+0.5),(float)(y_f+0.5)}};
  float corner_radii[]={(float)sqrt(pow(corners[0][0]-x_c,2)+
				    pow(corners[0][1]-y_c,2)),
			(float)sqrt(pow(corners[0][0]-x_c,2)+
				    pow(corners[1][1]-y_c,2)),
			(float)sqrt(pow(corners[1][0]-x_c,2)+
				    pow(corners[0][1]-y_c,2)),
			(float)sqrt(pow(corners[1][0]-x_c,2)+
				    pow(corners[1][1]-y_c,2))};

  float corner_phi[]={(float)atan2(corners[0][1]-y_c,
				   corners[0][0]-x_c),
		      (float)atan2(corners[1][1]-y_c,
				   corners[0][0]-x_c),
		      (float)atan2(corners[0][1]-y_c,
				   corners[1][0]-x_c),
		      (float)atan2(corners[1][1]-y_c,
				   corners[1][0]-x_c)};

  float deg=atan2(1,1)/45.0;
  float projected_x_coords[]={
    (float)(corner_radii[0]*cos(theta*deg-corner_phi[0])),
    (float)(corner_radii[1]*cos(theta*deg-corner_phi[1])),
    (float)(corner_radii[2]*cos(theta*deg-corner_phi[2])),
    (float)(corner_radii[3]*cos(theta*deg-corner_phi[3]))};
  float cs=cos(theta*deg);
  float sn=sin(theta*deg);

  if (g_verbose)
    cout << "corner_radii: " << 
      corner_radii[0] << " " << corner_radii[1] << " " <<
      corner_radii[2] << " " << corner_radii[3] << endl;

  if (g_verbose)
    cout << "proj_x_coords: " << 
      projected_x_coords[0] << " " << projected_x_coords[1] << " " <<
      projected_x_coords[2] << " " << projected_x_coords[3] << endl;
  
  float projected_x_lim[]={projected_x_coords[0],projected_x_coords[0]};
  for (int corner=1;corner<4;corner++) {
    projected_x_lim[0]=MIN(projected_x_lim[0],projected_x_coords[corner]);
    projected_x_lim[1]=MAX(projected_x_lim[1],projected_x_coords[corner]);
  }

  int n_bin=fp->m_info.px[roi_ix].nbin=
    (int)ceil(projected_x_lim[1]-projected_x_lim[0]);

  int *proj_x_hist=fp->m_info.px[roi_ix].proj_x_hist=
    (int*)malloc(n_bin*sizeof(int));
  int *proj_x_npix=fp->m_info.px[roi_ix].proj_x_npix=
    (int*)malloc(n_bin*sizeof(int));
  float *proj_x_coord=fp->m_info.px[roi_ix].proj_x_coord=
    (float*)malloc(n_bin*sizeof(float));
  float *proj_x_density=fp->m_info.px[roi_ix].proj_x_density=
    (float*)malloc(n_bin*sizeof(float));
      
  { // initialize.
    int b=n_bin;
    while (b--) proj_x_hist[b]=proj_x_npix[b]=0;
  }

  char *m=fp->m_info.map;
  //      int i,j;
  long tot=0L;

  for (int j=y_s;j<=y_f;j++) {
    for (int i=x_s;i<=x_f;i++) {
      int ind=j*fp->m_info.naxes[0]+i;
      int signal=m[ind];
      tot += signal;
      // figure out which histogram bin gets this contribution
      int bin_ix=floor(((i-x_c)*cs + (j-y_c)*sn)-projected_x_lim[0]);
      proj_x_hist[bin_ix] += signal;
      proj_x_npix[bin_ix] += 1;
    }
  }

  { // send out the histogram..
    char outfile[1024];
    sprintf(outfile,"%s_rix_%07d.roi",fp->filename,roi_ix);
    FILE *of=fopen(outfile,"w");
    for (int b=0;b<fp->m_info.px[roi_ix].nbin;b++) {
      proj_x_density[b]=(proj_x_hist[b]/((1.0)*proj_x_npix[b]));
      // should replace the x-axis (coord) with a better choice .. optionally,
      // since the increment is always unity, these can be turned into a
      // (float) starting point and a step size (if not unity)
      proj_x_coord[b]=(float)b;
      fprintf (of,"%d %d %d %g\n",
	       b,proj_x_npix[b],proj_x_hist[b],proj_x_density[b]);
    }
    fclose(of);
  }
  // dealocate here
  free(fp->m_info.px[roi_ix].proj_x_hist);
  free(fp->m_info.px[roi_ix].proj_x_npix);
  free(fp->m_info.px[roi_ix].proj_x_coord);
  free(fp->m_info.px[roi_ix].proj_x_density);

  fp->m_info.px[roi_ix].proj_x_hist=NULL;
  fp->m_info.px[roi_ix].proj_x_npix=NULL;
  fp->m_info.px[roi_ix].proj_x_coord=NULL;
  fp->m_info.px[roi_ix].proj_x_density=NULL;
}

void fits_report_error (int status) {
  fits_report_error(stderr,status);
  exit(1);
}
