'''
Created on Aug 22, 2012

@author: smirarab
'''
from multiprocessing import cpu_count, Pool, Lock
import copy
import traceback

class Job(object):
    '''
    The base class for a job that can be added to and executed form our JobPool.
    Subclasses of Job need to overwrite the 'run' method to perform their tasks. 
    A correct parallel way of running the job is by adding it to the JobPool.
     
    A Job object has  no direct access to its AsyncResult object, and thus has 
    no way of checking whether it has finished (for good reasons).
    AsyncResult object of a job, and more generally checking its status, should
    instead occur through the JobPool object. All job management tasks should
    occur through JobPool. job.result should be accessed only once one is sure t
    hat there the job has finished successfully (through JobPool).   

    A Job keeps a list of callback functions, which are called (in the order 
    they were added) when the job finishes running. Note that for addressing 
    some pickling issues, while a job is running, this list is copied into the
    JobPool object, and is set to None in the Job object. It should be noted
    that a Callback function is expected to never raise an exception. If an 
    exception is raised, it is simply printed out, but is otherwise totally
    ignored. To avoid deadlock, the status of the job is set to successful even
    if callbacks throw an exception. 
    
    In general, if the "run" method of a job modifies the state of the Job object, 
    these changes will NOT be propagated to other processes, including the parent
    process that runs the JobPool. The only thing that gets passed to the parent
    process is the "result" of a job. The result of the job is simply the return
    value of the "run" method. It gets automatically set to self.result once the
    job finishes running *successfully*. callback functions are also called only
    if Job finishes successfully, and therefore can assume that job.result is set.

    If "run" raises an exception, its job.result is never set.  
    
    Similarly, if the state of a Job object is modified after it is put in the 
    queue, the process that runs the Job might not see those changes.  
    
    Callback functions can access global state, but they should do so in a
    thread-safe fashion, because other job's callback functions could be trying 
    to access the same global state at the same time, giving rise to race 
    conditions. Use of a Lock, a Queue, or Pipe are recommended. 
    '''

    def __init__(self):
        '''
        Constructor
        '''
        self.callbacks = []
        self.add_call_Back(self._finished)
        self.result = None
        self.resultSet = False
        self.err = None

    def __call__(self):
        ''' This method makes this class a callable object'''
        return self.run()
        
    def run(self):
        '''
        This is the method that does the actual job. It should be overridden
        in subclasses. The return value is what gets set as self.result, and
        is passed from child process to the parent process (and should be
        therefore picklable)
        '''
        raise NotImplementedError;      
        
    def add_call_Back(self, function):
        ''' Adds a callback function to this job. The callback function should
        accept one parameter, which is going to be set to self.result.
        
        Callback functions should be thread-safe. If the access global state, 
        they should do with care. They should also return rather quickly, as
        they can block others.
        
        It is safe (to the best of our effort/knowledge) for a callback function
        to enqueue more jobs.
        '''
        if not hasattr(function, "__call__"):
            raise TypeError("Expecting a callable object. Instead getting a %s" %type(function))
        self.callbacks.append(function)
                    
    def _finished(self, results):
        '''
        This is a default callback function automatically added for all jobs. 
        It is the first function called when job finishes, and simply sets
        the result of the job.
        '''
        self.result = results
        self.resultSet = True
       

'''
JobPool is provided as a singleton class
'''
_jobPool = None
def JobPool(cpus = None):
    ''' The JobPool singelton instance is returned. 
    If the function is called for the first time, a new instance is created with
    the given number of cpus. Subsequent calls return the existing instances,
    if the number of requested cpus match the existing instance, or is None.
    If the number of requested cpus in subsequent calls is different from that
    of the original call, an error is raised. 
    '''
    global _jobPool
    if _jobPool is None: 
        _jobPool = _JobPool(cpus)
    if cpus is not None and cpus != _jobPool.cpus:
        raise Exception("An instance of Singleton class JobPool with %d cpus exists. You requested %d cpus" 
                        %( _jobPool.cpus, cpus)) 
    return _jobPool
    
class _JobPool:
    '''
    A singelton class, managing all Jobs. 
    Use JobPool() function to get the instance.
    Can be initialized only once, with a given number of cpus. 
    If number of cpus is left blank, all but two cpus (or only 1 cpu if less than three available) are used
    Jobs can be queued for execution. When a job finishes, its callback functions are called.
    One can choose to wait for all jobs in the queue to stop. 
    Scheduling and managing jobs
    '''
    
    def __init__(self, cpunum):
        self.cpus = cpunum if cpunum is not None else max([cpu_count()-2,1])
        self.pool = Pool(self.cpus)
        self._async_restults = {}
        self._callBack_lists = {}
        self.lock = Lock()
        
    def enqueue_job(self, job):

        if not isinstance(job, Job):
            raise TypeError("Expecting a Job object, instead getting a %s" %type(job))
        
        def call_back(job, result, callBackCopy):
            try:
                '''First restore callbacks attribute of the job'''
                job.callbacks = callBackCopy
                '''Then run all the callback in order, passing result as parameter'''
                for callback in callBackCopy:
                    apply(callback, [result])
            except Exception:
                # TODO: currently callback exceptions are simply ignored. 
                # What is a better solution?
                traceback.print_exc()
                
        ''' We need to backup callbacks of a Job. 
        The following line ensures that pickling issues do not arise'''
        callBackCopy = copy.copy(job.callbacks)
        job.callbacks = None
        
        ''' Add jobs to the pool, and save resulting AsyncResult function in a
        safe manner'''
        self.lock.acquire()
        result = self.pool.apply_async(job, callback = lambda res: call_back(job, res , callBackCopy))
        self._async_restults[job] = result
        self.lock.release()
        

    def wait_for_all_jobs(self, ignore_error=False):
        '''
        waits for all jobs in the queue to finish running. Makes attempt to 
        wait for newly added jobs as well.
        If ignore_error is false, this method stops after the first job with an
        error is discovered, and raises the error of that method. 
        Returns whether all finished jobs have finished successfully (in case
        ignore_error is set to True). 
        '''
        more = True
        while more:
            for result in self._async_restults.values():
                try:
                    result.get()
                except Exception as e:
                    if not ignore_error:
                        raise e
            ''' lock is acquired so that while we are iterating through the list
            no other job can add to the list (i.e. from their callbacks)'''            
            self.lock.acquire()
            more = False 
            hasFailed = False           
            for result in self._async_restults.values():
                if not result.ready():
                    more = True
                elif not result.successful():
                    if not ignore_error:
                        result.get()
                    hasFailed = True
            self.lock.release()
        
        return not hasFailed

    def is_job_running(self, job):
        '''Check to see if a given job has finished running'''
        return not self.get_job_asynch_result_object(job).ready()

    def get_job_asynch_result_object(self, job):
        ''' returns an instance of multiprocessing.pool.AsyncResult class,
        corresponding to the execution of the given job. This object can
        be used to check status of a job'''
        self.lock.acquire()
        if job not in self._async_restults.keys():
            raise KeyError("Job %s not in the queue yet." %str(job))
        resobj = self._async_restults[job]
        self.lock.release()
        return resobj
        
    def get_failed_jobs(self):
        ''' Returns a list of files that have failed (i.e. their run method has
        finished with an exception raised'''
        ret = []
        self.lock.acquire()
        for job, res in self._async_restults.iteritems():
            if res.ready() and not res.successful():
                ret.append(job)
        self.lock.release()
        return ret
    
    def get_job_error(self,job):
        '''
        If the given job has finished successfully, this method simply returns None.
        Otherwise, it raised the exception that terminated the given job 
        '''
        res = self.get_job_asynch_result_object(job)
        if not res.ready() or res.successful():
            return None
        else:
            try:
                res.get()
            except Exception as e:
                return e
            
    def get_all_job_errors(self):
        ''' Returns a list of all exceptions raised in all jobs executed in this
        pool. An empty list means no exceptions '''
        self.lock.acquire()
        jobs = self._async_restults.keys()
        self.lock.release()
        ret = [x for x in 
               [self.get_job_error(job) for job in jobs] 
                          if x is not None]
        return ret
        