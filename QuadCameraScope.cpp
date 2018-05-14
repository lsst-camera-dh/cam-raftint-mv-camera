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
#define _POSIX_C_SOURCE 200809L
#include <inttypes.h>
#include <stdio.h>
#include <time.h>

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
int g_save_roi_traces=0;
unsigned short **g_stacked_image;
void **g_normal_image;
unsigned short **g_pix_Xlut;
unsigned short **g_pix_Ylut;
int *g_stacked_frame_ix;
int g_filename_ringbuffer=0;
int *g_filename_ringbuffer_ix;
int g_savestack_nth_frame=0;
int g_save_rebin=0;
int g_save_mosaic=0;
// int g_device_map[]={1,0,3,2};
int g_device_map[]={0,2,1,3};
int g_device_flip_x[]={1,1,0,0};
int g_device_flip_y[]={0,0,1,1};
int g_pLock=0;
unsigned int *g_pMosaic=NULL;
int g_decimate=0;
int g_mosaic_rebin=4;
int g_mosaicNX=0;
int g_mosaicNY=0;
int g_mosaic_tile[]={2,2};
char **g_serials;
int g_bpp=8;
int g_ncameras=0;
long g_expotime=200000;
char g_lvl=0;
char g_signalhists=0;
time_t g_start_ts;
FILE *g_scopetrace_file;

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
  float proj_x_density_lims[2];
} proj_x;
  
typedef struct map_info_s {
  long   naxes[2];
  void   *map;
  int    n_px;
  proj_x **pxp;
} map_info;

typedef struct fits_s {
  char *filename;
  fitsfile *ff;
  map_info m_info;
  time_t   ts;
  struct fits_s *next_fits;
} fits;

void append_fits_request (fits **f,char *name,Request *thisreq,int camera_index,time_t ts);
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
  proj_x *px;
  struct roi_s *next_roi;
  // other items that are derived from the image in question and this extraction spec.
} roi;

void append_roi  (roi **r,char *roispec);
void report_rois (roi *r);

time_t msec_timestamp ();

static roi *g_roih=NULL;

float project_x (float x,float y,float cs, float sn,
		 float x_c,float y_c,float x_o);

float lvl_crossing_x (float norm_lvl,proj_x* px);

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

void adjust_exposure (unsigned int npix, void* pData, 
		      unsigned short dmask,
		      Device* pDev,time_t ts,int camera_index);


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
  switch (g_bpp) {
  case 8:
    ifc.pixelFormat.writeS("Mono8");
    break;
  case 10:
    if (strcmp("Mono10",ifc.pixelFormat.readS().c_str())!=0)
      ifc.pixelFormat.writeS("Mono10");
    break;
  case 12:
    if (strcmp("Mono12",ifc.pixelFormat.readS().c_str())!=0)
      ifc.pixelFormat.writeS("Mono12");
    break;
  default:
    break;
  }
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
       {"bpp",     required_argument, 0,             'B'},
       {"ring",    required_argument, 0,             'G'},
       {"stack",   required_argument, 0,             'S'},
       {"signalhist",no_argument,     0,             'g'},
       {"expoauto",required_argument, 0,             'a'},
       {"expotime",required_argument, 0,             'x'},
       {"save_raw",no_argument,       0,             's'},
       {"mosaic",  required_argument, 0,             'm'},
       {"roi",     required_argument, 0,             'r'},
       {"save_roi_traces",no_argument,0,             't'},
       {0,         0,                 0,             0}};

    int option_index=0;
    c=getopt_long(argc,argv,"tgbdsm:r:x:G:",
		  long_options,&option_index);

    if (c==-1) break;

    switch (c) {
    case 0:
      if (long_options[option_index].flag != 0) break;
      printf ("option %s",long_options[option_index].name);
      if (optarg) printf(" with arg %s",optarg);
      printf ("\n");
      break;

    case 'B':
      printf("option -B (bpp) with value `%s'\n",optarg);
      sscanf(optarg,"%d",&g_bpp);
      if ((g_bpp-8)*(g_bpp-10)*(g_bpp-12) != 0) {
	cout << "--bpp must be followed by a valid bits per pixel specification. choose from (8|10|12)." << endl;
	exit(1);
      }
      break;

    case 'a':
      printf("option -a (expoauto) with value `%s'\n",optarg);
      sscanf(optarg,"%c",&g_lvl);
      break;

    case 'x':
      printf("option -x (expotime) with value `%s'\n",optarg);
      sscanf(optarg,"%ld",&g_expotime);
      break;

    case 'S':
      printf("option -S (stack) with value `%s'\n",optarg);
      sscanf(optarg,"%d",&g_savestack_nth_frame);
      break;

    case 'b':
      printf("option -b (bin) with value `%s'\n",optarg);
      sscanf(optarg,"%d",&g_save_rebin);
      break;

    case 'G':
      printf("option -G (ringbuffer) with value `%s'\n",optarg);
      sscanf(optarg,"%d",&g_filename_ringbuffer);
      break;

    case 's':
      printf("option -s (save raw)\n");
      g_save_raw=1;
      break;

    case 'g':
      printf("option -g (signalhists)\n");
      g_signalhists=1;
      break;

    case 't':
      printf("option -t (save_roi_traces)\n");
      g_save_roi_traces=1;
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
  
  // get starting timestamp
  g_start_ts=msec_timestamp();

  // prepare scopetrace filename and open file, if g_roih isn't NULL
  if (g_roih!=NULL) {
    char outfile_scopetrace[1024];
    sprintf(outfile_scopetrace,"dat/scopetrace_%ld.qdp",g_start_ts);
    g_scopetrace_file=fopen(outfile_scopetrace,"w");
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

    if (g_serials==NULL) {
      cout << "can't allocate! gserials -- exiting." 
	   << endl; 
      exit(1);
    }
    g_normal_image=(void **)
      malloc(g_ncameras*sizeof(void*));
    if (g_normal_image==NULL) {
      cout << "can't allocate! g_normal_image -- exiting." 
	   << endl; 
      exit(1);
    }
    g_pix_Xlut=(unsigned short**)malloc(g_ncameras*sizeof(unsigned short*));
    g_pix_Ylut=(unsigned short**)malloc(g_ncameras*sizeof(unsigned short*));
    if ((g_pix_Xlut==NULL) || (g_pix_Ylut==NULL)) {
      cout << "can't allocate! g_pix_*lut -- exiting." 
	   << endl; 
      exit(1);
    }
    
    if (g_filename_ringbuffer != 0) {
      g_filename_ringbuffer_ix=(int*)
	malloc(g_ncameras*sizeof(int));
      for (int j=0;j<g_ncameras;j++)
	g_filename_ringbuffer_ix[j]=0;
      // ready to use.
    }

    if (g_savestack_nth_frame != 0) {
      g_stacked_image=(unsigned short**)
	malloc(g_ncameras*sizeof(unsigned short*));
      if (g_stacked_image==NULL) {
	cout << "can't allocate! g_stacked_image -- exiting." 
	     << endl; 
	exit(1);
      }
      g_stacked_frame_ix=(int*)malloc(g_ncameras*sizeof(int));
    }

    for (int j=0;j<g_ncameras;j++) {
      g_normal_image[j]=NULL;
      g_pix_Xlut[j]=g_pix_Ylut[j]=NULL;
      g_serials[j]=(char*)malloc(1024*sizeof(char));
      if (g_serials[j]==NULL) {
	cout << "can't allocate! gserials[j] -- exiting." << endl;
	exit(1);
      }
      if (g_savestack_nth_frame) {
	g_stacked_image[j]=NULL; // g_stacked_image allocated below
	g_stacked_frame_ix[j]=0;
      }
    }
    g_mosaic_tile[1] = ceil(g_ncameras / g_mosaic_tile[0]);
    for( unsigned int i = 0; i < devCnt; i++ )	// populate globals for use later
      sprintf(g_serials[i],devMgr[i]->serial.read().c_str());
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
  if (! thisreq->isOK()) return;
  //  unsigned int ts = thisreq->infoTimeStamp_us.read();
  time_t ts = msec_timestamp();

  unsigned int xt = thisreq->infoExposeTime_us.read();
  string ser = pDev->serial.read();
  //    unsigned int frame_ix = thisreq->infoFrameID.read();
  // storage type of pData depends on g_bpp setting
  void* pData = thisreq->imageData.read();
  unsigned int iWidth = thisreq->imageWidthTotal.read();
  unsigned int iHeight = thisreq->imageHeightTotal.read();
  unsigned int npix=iWidth*iHeight;

  int camera_index=g_ncameras;
  while (camera_index-- && 
	 strcmp(g_serials[camera_index],pDev->serial.read().c_str()));
  // and use the mapping..
  camera_index=g_device_map[camera_index];


  // choose the right dmask..
  unsigned short dmask;
  switch (g_bpp) {
  case 8:
    dmask=0x00FF;
    break;
  case 10:
    dmask=0x03FF;
    break;
  case 12:
    dmask=0x0FFF;
    break;
  default:
    cout << "g_bpp is " << g_bpp << ". this shouldn't happen." << endl;
    exit(1);
    break;
  }

  // allocate the normal buffer if not already done..
  if (g_normal_image[camera_index]==NULL) {
    // then allocate:
    g_normal_image[camera_index]=
      (void*)malloc(iWidth*iHeight*
		    ((g_bpp==8)?sizeof(char):sizeof(unsigned short)));
    
    if (g_normal_image[camera_index]==NULL) {
      cout << "can't allocate! g_normal_image[camera_index] -- exiting." << endl; 
      exit(1);
    }

    g_pix_Xlut[camera_index]=(unsigned short*)malloc(iWidth*sizeof(unsigned short));
    g_pix_Ylut[camera_index]=(unsigned short*)malloc(iHeight*sizeof(unsigned short));
    if ((g_pix_Xlut[camera_index]==NULL) || (g_pix_Ylut[camera_index]==NULL)) {
      cout << "can't allocate! g_pix_*lut[camera_index] -- exiting." << endl; 
      exit(1);
    }
    // prepare reference array (map) to correct mapping for
    // ds9 display etc. (using g_device_flip_x[] and g_device_flip_y[].
    // here fill in the g_pix_addr_map for mapping incoming pixels.
    int i;
    i=iWidth;
    while (i--)
      g_pix_Xlut[camera_index][(g_device_flip_x[camera_index])?iWidth-1-i:i]=i;
    i=iHeight;
    while (i--)
      g_pix_Ylut[camera_index][(g_device_flip_y[camera_index])?iHeight-1-i:i]=i;

  } // won't deallocate this buffer or LUTs on subroutine exit, will use throughout program


  // map the incoming pData onto g_normal_image and then reassign pData to g_normal_image
  {
    void *d=g_normal_image[camera_index];
    unsigned short *xlut=g_pix_Xlut[camera_index];
    unsigned short *ylut=g_pix_Ylut[camera_index];
    int i=iWidth*iHeight;
    switch (g_bpp) {
    case 8:
      while (i--) 
	((char*)d)[ylut[i/iWidth]*iWidth+xlut[i%iWidth]] = 
	  (dmask & ((char*)pData)[i]);
      break;
    case 10:
    case 12:
      while (i--) 
	((unsigned short*)d)[ylut[i/iWidth]*iWidth+xlut[i%iWidth]] = 
	  (dmask & ((unsigned short*)pData)[i]);
      break;
    default:
      break;
    }
    // now set pData to g_normal_image:
    pData=g_normal_image[camera_index];
  }

  char outfile_raw[1024],outfile_raw_stage[1024],outfile_raw_save[1024];

  if (g_filename_ringbuffer != 0) {
    sprintf(outfile_raw,"dat/BF%1d_rng%04d_%u.fits",camera_index,(g_filename_ringbuffer_ix[camera_index]%g_filename_ringbuffer),xt);
    sprintf(outfile_raw_stage,"dat/BF%1d_rng%04d_%u.fits.partial",camera_index,(g_filename_ringbuffer_ix[camera_index]%g_filename_ringbuffer),xt);
    sprintf(outfile_raw_save,
	    "dat/BF%1d_rng%04d_%u.fits",
	    camera_index,(g_filename_ringbuffer_ix[camera_index]%g_filename_ringbuffer),xt);
  } else {
    // sprintf(outfile_raw,
    //  "dat/BF_%s_%ld_%u.fits",ser.c_str(),ts,xt);
    sprintf(outfile_raw,
	    "dat/BF%1d_%ld_%u.fits",camera_index,ts,xt);
    sprintf(outfile_raw_stage,"dat/BF%1d_%ld_%u.fits.partial",camera_index,ts,xt);
    sprintf(outfile_raw_save,"dat/BF%1d_%ld_%u.fits",camera_index,ts,xt);
  }

  if (g_save_raw) {
    int status=0;
    int stage_rename=1;
    fitsfile *ff;
    // fills disk space quickly and can maintain ~ 33 fps * 5Mbyte on 4 cameras simultaneously without impacting frame rate
    if (g_verbose) cout << " (" << outfile_raw << ") ";
    if (stage_rename) {
      fits_create_file(&ff,outfile_raw_stage,&status);
    } else {
      fits_create_file(&ff,outfile_raw,&status);
    }
    fits_report_error(status);
    long naxis=2;
    long naxes[]={(long)iWidth,(long)iHeight};
    switch (g_bpp) {
    case 8:
      fits_create_img(ff,BYTE_IMG,naxis,naxes,&status);
      break;
    case 10:
    case 12:
      fits_create_img(ff,USHORT_IMG,naxis,naxes,&status);
      break;
    default:
      cout << "unexpected bpp.. this shouldn't happen!" << endl;
      exit(1);
      break;
    }
    fits_report_error(status);
    fits_write_img(ff,((g_bpp==8)?TBYTE:TUSHORT),
		   1,naxes[0]*naxes[1],
		   (void*)g_normal_image[camera_index],&status);
    fits_report_error(status);
    fits_close_file(ff,&status);
    fits_report_error(status);
    if (stage_rename) {
      int rc;
      remove(outfile_raw_save); // might fail if file doesn't exist.
      rc=rename(outfile_raw_stage,outfile_raw_save);
      if (rc) {
	perror("error renaming file.\n"); 
	exit(1);
      }
    }
  }
  
  if (g_save_rebin!=0) { 
    // save another representation of the image, 
    // this time decimated by g_save_rebin x g_save_rebin for smaller size

    char outfile_rebin[1024];
    void *pData_rebin;
    int status=0;
    fitsfile *ff;
    long naxis=2;
    long naxes[]={(long)floor(iWidth/g_save_rebin),(long)floor(iHeight/g_save_rebin)};
      
    if (g_decimate) {
      if (g_filename_ringbuffer != 0) {
	sprintf(outfile_rebin,
		"!dat/BF%1d_dm%02d_rng%04d_%u.fits",camera_index,
		g_save_rebin,
		g_filename_ringbuffer_ix[camera_index]%
		g_filename_ringbuffer,xt);
      } else {
	//	  sprintf(outfile_rebin,"dat/BF_%s_dm%02d_%u_%u.fits",
	//		  ser.c_str(),g_save_rebin,ts,xt);
	sprintf(outfile_rebin,"dat/BF%1d_dm%02d_%ld_%u.fits",
		camera_index,g_save_rebin,ts,xt);
      }
    } else {
      if (g_filename_ringbuffer != 0) {
	sprintf(outfile_rebin,
		"!dat/BF%1d_rb%02d_rng%04d_%u.fits",camera_index,
		g_save_rebin,
		g_filename_ringbuffer_ix[camera_index]%
		g_filename_ringbuffer,xt);
      } else {
	sprintf(outfile_rebin,"dat/BF_%s_rb%02d_%ld_%u.fits",
		ser.c_str(),g_save_rebin,ts,xt);
	sprintf(outfile_rebin,"dat/BF%1d_rb%02d_%ld_%u.fits",camera_index,
		g_save_rebin,ts,xt);
      }
    }
    // pData_rebin type no longer depends on whether raw data will be decimated or rebinned.
    pData_rebin=(void*)malloc(naxes[0]*naxes[1]*sizeof(unsigned short));

    if (pData_rebin==NULL) {
      cout << "can't allocate! pData_rebin -- exiting." << endl;
      exit(1);
    }

    if (g_verbose)
      cout << " (" << outfile_rebin << ") " << endl;

    fits_create_file(&ff,outfile_rebin,&status);
    fits_report_error(status);

    if (g_decimate) {
      // decimate, not rebin
      int ix=naxes[0]*naxes[1];
      switch (g_bpp) {
      case 8:
	while (ix--) {
	  int i=(ix%naxes[0])*g_save_rebin;
	  int j=(ix/naxes[0])*g_save_rebin;
	  ((unsigned short*)pData_rebin)[ix]=
	    (unsigned short)((char*)(pData))[j*iWidth+i];
	}
	break;
      case 10:
      case 12:
	while (ix--) {
	  int i=(ix%naxes[0])*g_save_rebin;
	  int j=(ix/naxes[0])*g_save_rebin;
	  ((unsigned short*)pData_rebin)[ix]=
	    ((unsigned short*)(pData))[j*iWidth+i];
	}
	break;
      default:
	break;
      }

    } else {
      // rebin, not decimate
      int ix=naxes[0]*naxes[1];
      while (ix--) ((unsigned short*)pData_rebin)[ix]=0;
      ix=iWidth*iHeight;
      switch (g_bpp) {
      case 8:
	while (ix--) {
	  int i=(ix%iWidth)/g_save_rebin;
	  int j=(ix/iWidth)/g_save_rebin;
	  ((unsigned short*)pData_rebin)[j*naxes[0]+i] += 
	    (unsigned short)(dmask & ((char*)(pData))[ix]);
	}
	break;
      case 10:
      case 12:
	while (ix--) {
	  int i=(ix%iWidth)/g_save_rebin;
	  int j=(ix/iWidth)/g_save_rebin;
	  ((unsigned short*)pData_rebin)[j*naxes[0]+i] +=
	    (dmask & ((unsigned short*)(pData))[ix]);
	}
	break;
      default:
	break;
      }
    }

    fits_create_img(ff,USHORT_IMG,naxis,naxes,&status);
    fits_report_error(status);
    fits_write_img(ff,TUSHORT,1,naxes[0]*naxes[1],
		   (void*)pData_rebin,&status);
    fits_report_error(status);
    fits_close_file(ff,&status);
    fits_report_error(status);

    if (g_verbose)
      cout << "about to free pData_rebin " << (void*)pData_rebin << endl;
    free(pData_rebin);
    pData_rebin=NULL;
  }

  if (g_savestack_nth_frame != 0) {

    { // initialize if necessary
      if (g_stacked_image[camera_index]==NULL) {
	g_stacked_image[camera_index]=(unsigned short*)malloc(iWidth*iHeight*sizeof(unsigned short));
	if (g_stacked_image[camera_index]==NULL) {
	  cout << "can't allocate! g_stacked_image[camera_index] -- exiting." << endl; 
	  exit(1);
	}
	int i=iWidth*iHeight;
	// initialize
	while (i--) g_stacked_image[camera_index][i]=0; 
      }
    }

    int i=iWidth*iHeight;
    switch (g_bpp) {
    case 8:
      while (i--)
	g_stacked_image[camera_index][i]+=
	  (unsigned short)(dmask & ((char*)pData)[i]);
      break;
    case 10:
    case 12:
      while (i--)
	g_stacked_image[camera_index][i]+=
	  (dmask & ((unsigned short*)pData)[i]);
      break;
    }
    
    g_stacked_frame_ix[camera_index]++;
    
    if (g_stacked_frame_ix[camera_index] % g_savestack_nth_frame == 0) {
      char outfile_stack[1024];
      sprintf(outfile_stack,"dat/BF_%s_stack%02d_%ld_%u.fits",
	      ser.c_str(),g_savestack_nth_frame,ts,xt);
      sprintf(outfile_stack,"dat/BF%1d_stack%02d_%ld_%u.fits",camera_index,
	      g_savestack_nth_frame,ts,xt);
      if (g_filename_ringbuffer != 0) {
	sprintf(outfile_stack,"!dat/BF%1d_stack%02d_rng%04d_%u.fits",camera_index,
		g_savestack_nth_frame,
		(g_stacked_frame_ix[camera_index]/g_savestack_nth_frame-1)%g_filename_ringbuffer,
		xt);
      }
      // save the image and zero it out again.
      int status=0;
      fitsfile *ff;
      if (g_verbose)
	cout << " (" << outfile_stack << ") " << endl;
      fits_create_file(&ff,outfile_stack,&status);
      fits_report_error(status);
      long naxis=2;
      long naxes[]={(long)iWidth,(long)iHeight};
      fits_create_img(ff,USHORT_IMG,naxis,naxes,&status);
      fits_report_error(status);
      fits_write_img(ff,TUSHORT,1,naxes[0]*naxes[1],
		     (void*)g_stacked_image[camera_index],
		     &status);
      fits_report_error(status);
      fits_close_file(ff,&status);
      fits_report_error(status);
      i=iWidth*iHeight;
      while (i--) g_stacked_image[camera_index][i] = 0;
    }
  }

  if (g_save_mosaic) {
    // figure out which sensor this is (by serial) to determine its place
    if (g_verbose)
      cout << "got camera index: " << camera_index << endl;
    // initialize if necessary, using g_pLock etc.
    while (g_pLock) {
      // do nothing until lock is released
      //	if (g_verbose)
      //	  cout << " waiting .. (camera_index=" << camera_index << ")";
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
      if (g_pMosaic==NULL) {cout << "can't allocate! g_pMosaic -- exiting." << endl; 
	exit(1);}
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
	//	  if (g_verbose)
	//	    cout << " waiting .. (camera_index=" << camera_index << ")";
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
      switch (g_bpp) {
      case 8:
	while (i--) {
	  int x=startx+(i%iWidth)/g_mosaic_rebin;
	  int y=starty+(i/iWidth)/g_mosaic_rebin;
	  g_pMosaic[y*g_mosaicNX+x] += 
	    (unsigned short)(dmask & ((char*)(pData))[i]);
	}
	break;
      case 10:
      case 12:
	while (i--) {
	  int x=startx+(i%iWidth)/g_mosaic_rebin;
	  int y=starty+(i/iWidth)/g_mosaic_rebin;
	  g_pMosaic[y*g_mosaicNX+x] += 
	    (dmask & ((unsigned short*)(pData))[i]);
	}
	break;
      default:
	break;
      }

      // save the mosaic if appropriate (do this only for a specific camera index case
      if (camera_index==0) {
	// make a copy of the data and unlock for other camera threads to proceed
	unsigned int* mosaic_copy=(unsigned int*)malloc(g_mosaicNX*g_mosaicNY*sizeof(unsigned int));
	if (mosaic_copy==NULL) {
	  cout <<"can't allocate! mosaic_copy -- exiting."<<endl; 
	  exit(1);
	}
	memcpy(mosaic_copy,g_pMosaic,g_mosaicNX*g_mosaicNY*sizeof(unsigned int));
	g_pLock=0;
	// unfinished business of saving the fits file and freeing up the copy image.
	int status=0;
	char outfile_mosaic[1024];
	fitsfile *ff;
	if (g_filename_ringbuffer != 0) {
	  sprintf(outfile_mosaic,
		  "!dat/BF_%s_rng%04d_%u.fits",
		  "mosaic",
		  (g_filename_ringbuffer_ix[camera_index]%g_filename_ringbuffer),
		  xt);
	} else {
	  sprintf(outfile_mosaic,"dat/BF_%s_%ld_%u.fits","mosaic",ts,xt);
	}
	if (g_verbose)
	  cout << " (" << outfile_mosaic << ") " << endl;
	  
	fits_create_file(&ff,outfile_mosaic,&status);
	fits_report_error(status);
	long naxis=2;
	long naxes[]={(long)g_mosaicNX,(long)g_mosaicNY};
	fits_create_img(ff,USHORT_IMG,naxis,naxes,&status);
	fits_report_error(status);
	fits_write_img(ff,TUINT,1,naxes[0]*naxes[1],(void*)mosaic_copy,&status);
	fits_report_error(status);
	fits_close_file(ff,&status);
	fits_report_error(status);
	free(mosaic_copy);
	mosaic_copy=NULL;
      } else {
	// if other camera_index, done with contributing to g_pMosaic
	g_pLock=0;
      }
    }
  }

  if (0) {
    cout << "timestamp is: " << ts << " camera_index is: " << camera_index << endl;
    cout << "dmask: " << dmask << " g_normal_image[camera_index]: " << g_normal_image[camera_index] << endl;
    cout << "g_pix_Xlut[camera_index] " << g_pix_Xlut[camera_index] << endl;
    cout << "g_pix_Ylut[camera_index] " << g_pix_Ylut[camera_index] << endl;
  }
  
  if (g_lvl || g_signalhists) 
    adjust_exposure(npix,pData,dmask,pDev,ts,camera_index);

  // finally compute rois on pData
  if (g_roih != NULL) {
    // append_fits_request gets its information & buffer from memory, not from the 
    // output file just saved (if in fact it was saved). name of fits file (outfile)
    // is used primarily for display purposes.
    // "fitsh" a single element fits list (be sure to deallocate after use)
    int camera_index=g_ncameras;
    while (camera_index-- && strcmp(g_serials[camera_index],pDev->serial.read().c_str()));
    camera_index=g_device_map[camera_index];
    fits *fitsh=NULL;
    append_fits_request(&fitsh,outfile_raw,thisreq,camera_index,ts); 
    compute_roi_fits_request(g_roih,fitsh,camera_index);
    destroy_fits_request(&fitsh);
  }

  if (g_filename_ringbuffer != 0)
    g_filename_ringbuffer_ix[camera_index]++;
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
    fp->m_info.pxp=(proj_x**)malloc(fp->m_info.n_px*sizeof(proj_x*));
    if (fp->m_info.pxp==NULL) {
      cout << "can't allocate! fp->m_info.pxp -- exiting." << endl;
      exit(1);
    }

    rp=r;
    roi_ix=0;
    while (rp!=NULL) {
      // important - apply the same filter as above here
      int thiscamera=rp->camera_index;
      if (cix == thiscamera) {
	fp->m_info.pxp[roi_ix]=rp->px;
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
	// do nothing, already deallocated
      }
      free(fp->m_info.pxp);
      fp->m_info.pxp=NULL;
    }
    fp=fp->next_fits;
  }
}

void append_roi (roi **r,char *roispec) {
  roi *new_roi=(roi*)malloc(sizeof(roi));
  if (new_roi==NULL) {
    cout << "can't allocate! new_roi -- exiting." << endl; 
    exit(1);
  }

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

  // allocate & initialize proj_x structure here.
  new_roi->px=(proj_x*)malloc(sizeof(proj_x));
  new_roi->px->nbin=0;
  new_roi->px->proj_x_hist=NULL;
  new_roi->px->proj_x_npix=NULL;
  new_roi->px->proj_x_coord=NULL;
  new_roi->px->proj_x_density=NULL;

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
    fp=NULL;
  }
}

void append_fits_request (fits **f,char *name,Request *thisreq,int camera_index,time_t ts) {
  fits *new_fits=(fits*)malloc(sizeof(fits));
  if (new_fits==NULL) {
    cout << "can't allocate! new_fits -- exiting." << endl; 
    exit(1);
  }
  new_fits->filename=name;
  new_fits->ts=ts;
  // copy the map over from the request (already in memory)
  new_fits->m_info.naxes[0]=thisreq->imageWidthTotal.read();
  new_fits->m_info.naxes[1]=thisreq->imageHeightTotal.read();
  // don't do this, the pixels may be in the wrong order
  // new_fits->m_info.map=(char*)thisreq->imageData.read(); // isn't malloc, don't free
  // do this:
  new_fits->m_info.map=g_normal_image[camera_index];
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
  int   cix=rp->camera_index;
  // avoid boundaries that violate limits
  int x_s = (int)((x0<0)?0:x0);
  int x_f = (int)((x1<fp->m_info.naxes[0]-1)?x1:fp->m_info.naxes[0]-1);
  int y_s = (int)((y0<0)?0:y0);
  int y_f = (int)((y1<fp->m_info.naxes[1]-1)?y1:fp->m_info.naxes[1]-1);
  // center coordinates for this roi
  float x_c,y_c,x_offset;
  // actually use the intended center, that's where this needs to rotate about

  if (0) { // rotation point about center of requested ROI (which might not be the same as center of available data)
    // original way:
    // x_c=(x_s+x_f)/2.0;
    // y_c=(y_s+y_f)/2.0;
    x_c = (x0+x1)/2.0;
    y_c = (y0+y1)/2.0;
    x_offset = x_c;
  } else { // rotation point about center of field
    x_c=(0+fp->m_info.naxes[0])/2.0;
    y_c=(0+fp->m_info.naxes[1])/2.0;
    x_offset = x_c;
  }

  // find the limits (in pixels) when the roi is rotated as prescribed
  float corners[][2]={{(float)(x_s-0.5),(float)(y_s-0.5)},
		      {(float)(x_f+0.5),(float)(y_f+0.5)}};

  float deg=atan2(1,1)/45.0;
  float cs=cos(theta*deg);
  float sn=sin(theta*deg);

  float proj_corners[]={
    project_x(corners[0][0],corners[0][1],cs,sn,x_c,y_c,x_offset),
    project_x(corners[0][0],corners[1][1],cs,sn,x_c,y_c,x_offset),
    project_x(corners[1][0],corners[0][1],cs,sn,x_c,y_c,x_offset),
    project_x(corners[1][0],corners[1][1],cs,sn,x_c,y_c,x_offset)
  };

  if (g_verbose)
    cout << "proj_x_corners: " << proj_corners[0] << " " << proj_corners[1] << " " << proj_corners[2] << " " << proj_corners[3] << endl;
  
  float projected_x_lim[]={proj_corners[0],proj_corners[0]};
  for (int corner=1;corner<4;corner++) {
    projected_x_lim[0]=MIN(projected_x_lim[0],proj_corners[corner]);
    projected_x_lim[1]=MAX(projected_x_lim[1],proj_corners[corner]);
  }

  int n_bin=fp->m_info.pxp[roi_ix]->nbin=
    (int)floor(projected_x_lim[1]-projected_x_lim[0]);
  int *proj_x_hist=fp->m_info.pxp[roi_ix]->proj_x_hist;
  if (proj_x_hist==NULL)
    proj_x_hist=fp->m_info.pxp[roi_ix]->proj_x_hist=
      (int*)malloc(n_bin*sizeof(int));
  if (proj_x_hist==NULL) {
    cout << "can't allocate! proj_x_hist -- exiting." << endl; 
    exit(1);
  }

  int *proj_x_npix=fp->m_info.pxp[roi_ix]->proj_x_npix;
  if (proj_x_npix==NULL)
    proj_x_npix=fp->m_info.pxp[roi_ix]->proj_x_npix=
      (int*)malloc(n_bin*sizeof(int));
  if (proj_x_npix==NULL) {
    cout << "can't allocate! proj_x_npix -- exiting." << endl; 
    exit(1);
  }

  float *proj_x_coord=fp->m_info.pxp[roi_ix]->proj_x_coord;
  if (proj_x_coord==NULL)
    proj_x_coord=fp->m_info.pxp[roi_ix]->proj_x_coord=
      (float*)malloc(n_bin*sizeof(float));
  if (proj_x_coord==NULL) {
    cout << "can't allocate! proj_x_coord -- exiting." << endl; 
    exit(1);
  }

  float *proj_x_density=fp->m_info.pxp[roi_ix]->proj_x_density;
  if (proj_x_density==NULL)
    proj_x_density=fp->m_info.pxp[roi_ix]->proj_x_density=
      (float*)malloc(n_bin*sizeof(float));
  if (proj_x_density==NULL) {
    cout << "can't allocate! proj_x_density -- exiting." << endl;
    exit(1);
  }


  { // initialize.
    int b=n_bin;
    while (b--) proj_x_hist[b]=proj_x_npix[b]=0;
  }

  void *m=fp->m_info.map;
  //      int i,j;
  long tot=0L;


  for (int j=y_s;j<=y_f;j++) {
    for (int i=x_s;i<=x_f;i++) {
      int ind=j*fp->m_info.naxes[0]+i;
      int signal=((g_bpp==8)?(int)((char*)m)[ind]:(int)((unsigned short*)m)[ind]);
      tot += signal;
      // figure out which histogram bin gets this contribution
      int bin_ix=floor(((i-x_c)*cs + (j-y_c)*sn)-projected_x_lim[0]);
      bin_ix=floor(project_x((float)i,(float)j,cs,sn,x_c,y_c,x_offset)
		   -projected_x_lim[0]);
      proj_x_hist[bin_ix] += signal;
      proj_x_npix[bin_ix] += 1;
    }
  }

  // fill out the proj_x_density and proj_x_coord arrays. parasitically record min/max density
  // for estimating edge sharpness.
  // point density_lims to structure address
  float *proj_x_density_lims=
    fp->m_info.pxp[roi_ix]->proj_x_density_lims;

  for (int b=0;b<fp->m_info.pxp[roi_ix]->nbin;b++) {
    proj_x_density[b]=(proj_x_hist[b]/((1.0)*proj_x_npix[b]));
    proj_x_coord[b] = projected_x_lim[0] + b;
    if (b==0) {
      proj_x_density_lims[0]=proj_x_density_lims[1]=proj_x_density[b];
    } else {
      proj_x_density_lims[0]=MIN(proj_x_density_lims[0],proj_x_density[b]);
      proj_x_density_lims[1]=MAX(proj_x_density_lims[1],proj_x_density[b]);
    }
  }

  if (1) {
    // compute projected x positions for the following density levels: (can use lvl_crossing_x for angle finding and
    // for time series/scope results)
    float xval;
    int nlvls=3;
    float lvls[]={0.1,0.5,0.9};
    char scope_trace[2048];
    sprintf(scope_trace,"%ld %d %d %f ",fp->ts-g_start_ts,cix,roi_ix,theta);
    
    for (int lvl_ix=0;lvl_ix<nlvls;lvl_ix++) {
      xval=lvl_crossing_x(lvls[lvl_ix],fp->m_info.pxp[roi_ix]);
      cout << "crossing @ " << lvls[lvl_ix] << " x = " << xval << endl;
      sprintf(scope_trace,"%s %f ",scope_trace,xval);
    }
    fprintf(g_scopetrace_file,"%s\n",scope_trace);
    fflush(g_scopetrace_file);
  }
  

  if (g_save_roi_traces) { // send out the histogram..
    char outfile[1024];
    sprintf(outfile,"%s_rix_%07d.roi",fp->filename,roi_ix);
    cout << "target filename for rix: " << outfile << endl;
    FILE *of=fopen(outfile,"w");
    for (int b=0;b<fp->m_info.pxp[roi_ix]->nbin;b++) 
      fprintf (of,"%g %g\n",proj_x_coord[b],proj_x_density[b]);
    fclose(of);
  }
}

void fits_report_error (int status) {
  if (status) {
    fits_report_error(stderr,status);
    exit(1);
  }
  return;
}

float project_x (float x,float y,float cs, float sn,
		 float x_c,float y_c,float x_o) {
  // seems that in all cases x_o would be equal to x_c 
  // but leave this as a separate parameter for now..
  return((x-x_c)*cs+(y-y_c)*sn+x_o);
}

float lvl_crossing_x (float norm_lvl,proj_x* px) {
  float *proj_x_density=px->proj_x_density;
  float *proj_x_coord=px->proj_x_coord;
  float *proj_x_density_lims=px->proj_x_density_lims;
  int nbin=px->nbin;
  float lvl=proj_x_density_lims[0]+norm_lvl*
    (proj_x_density_lims[1]-proj_x_density_lims[0]);

  int bin_ix=1;
  while ((bin_ix < nbin) && 
	 ((proj_x_density[bin_ix]-lvl)*
	  (lvl-proj_x_density[bin_ix-1])<0)) 
    bin_ix++;

  if (bin_ix>nbin-1) bin_ix=nbin-1;
  
  float lowerlim=proj_x_coord[bin_ix-1];
  float partial=(lvl-proj_x_density[bin_ix-1])/
    (proj_x_density[bin_ix]-proj_x_density[bin_ix-1]);
  float step=(proj_x_coord[bin_ix]-proj_x_coord[bin_ix-1]);
  return(lowerlim+partial*step);
}

time_t msec_timestamp () {
  // get the timestamp
  long ms; // milliseconds
  time_t s; //seconds
  struct timespec spec;

  clock_gettime(CLOCK_REALTIME,&spec);
  s = spec.tv_sec;
  ms=round(spec.tv_nsec/1.0e6); // convert nanoseconds to milliseconds
  if (ms>999) {
    s++;
    ms = 0;
  }
  //  cout << "Current time: " << (intmax_t)s << " ms " << ms << endl;
  return(s*1000+ms);
}

void adjust_exposure (unsigned int npix,void* pData,
		      unsigned short dmask,
		      Device *pDev,time_t ts,int camera_index) {

  if (!g_lvl && !g_signalhists) return;

  GenICam::AcquisitionControl ac( pDev );
  string ser = pDev->serial.read();
  unsigned long i;
  unsigned long tot_pix=0L,tot_sig=0L;

  int hist[dmask],chist[dmask];
  i=dmask;
  while (i--) hist[i]=0;

  for (i=0L;i<npix;i++) {
    unsigned short dat=
      ((g_bpp==8)?
       (unsigned short)((char*)(pData))[i]:
       ((unsigned short*)(pData))[i]);
    hist[dat]++;
    tot_sig += dat;
    tot_pix++;
  }

  // float avg_sig = tot_sig/((float)(tot_pix));
  // get quantile levels for adjusting exposure times

  const int nlvls=3;
  int ilvl;
  float lvls[]={0.1,0.5,0.9};
  float sig_lvls[nlvls];

  for (ilvl=0;ilvl<nlvls;ilvl++) {
    sig_lvls[ilvl]=-1.0;
  }

  for (i=0L;i<dmask;i++) {
    chist[i]=(i==0)?hist[i]:chist[i-1]+hist[i];
  }

  for (ilvl=0;ilvl<nlvls;ilvl++) {
    float lvl=lvls[ilvl]*npix;
    sig_lvls[ilvl]=0;
    for (i=1L;i<dmask;i++) {
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
  if (g_signalhists) { // and save the histogram 

    char histfile[1024];
    sprintf(histfile,"dat/BF_%s_%ld.hist",ser.c_str(),ts);
    sprintf(histfile,"dat/BF%1d_%ld.hist",camera_index,ts);
    FILE *fp;
    fp=fopen(histfile,"w");
    for (i=1L;i<dmask;i++) {
      fprintf(fp,"%d %g\n",(int)(i),chist[i]/(float)npix);
    }      
    fclose(fp);
  }
}

