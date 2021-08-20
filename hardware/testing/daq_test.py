import PyDAQmx as daq
import numpy as np

task_handle = daq.TaskHandle()
samps = 3000
test_array = np.full((samps,), -1, dtype=np.uint32)

daq.DAQmxCreateTask('count_test', daq.byref(task_handle))

daq.DAQmxCreateCICountEdgesChan(task_handle,
                                '/Dev1/Ctr0',
                                'Counter test',
                                daq.DAQmx_Val_Rising,
                                0,
                                daq.DAQmx_Val_CountUp)


daq.DAQmxCfgSampClkTiming(task_handle,
                          '100kHzTimebase',
                          100000,
                          daq.DAQmx_Val_Rising,
                          daq.DAQmx_Val_FiniteSamps,
                          samps
                          )


daq.DAQmxSetCICountEdgesTerm(task_handle,
                             '/Dev1/Ctr0',
                             '/Dev1/PFI0')

done = daq.bool32()
daq.DAQmxStartTask(task_handle)
daq.DAQmxGetTaskComplete(task_handle, done)

n_read_samples = daq.int32()
daq.DAQmxReadCounterU32(task_handle,
                        samps,
                        500.,
                        test_array,
                        len(test_array),
                        daq.byref(n_read_samples),
                        None
                        )

done = daq.bool32()
daq.DAQmxGetTaskComplete(task_handle, done)

daq.DAQmxStopTask(task_handle)

daq.DAQmxClearTask(task_handle)

print(test_array[-1] * 100000/samps)