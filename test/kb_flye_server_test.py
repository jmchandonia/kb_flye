# -*- coding: utf-8 -*-
import os
import time
import unittest
from configparser import ConfigParser

import shutil
import requests

from kb_flye.kb_flyeImpl import kb_flye
from kb_flye.kb_flyeServer import MethodContext
from kb_flye.authclient import KBaseAuth as _KBaseAuth
from installed_clients.ReadsUtilsClient import ReadsUtils
from installed_clients.AssemblyUtilClient import AssemblyUtil
from installed_clients.WorkspaceClient import Workspace
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.AbstractHandleClient import AbstractHandle as HandleService

class kb_flyeTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.token = os.environ.get('KB_AUTH_TOKEN', None)
        config_file = os.environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('kb_flye'):
            cls.cfg[nameval[0]] = nameval[1]
        # Getting username from Auth profile for token
        authServiceUrl = cls.cfg['auth-service-url']
        auth_client = _KBaseAuth(authServiceUrl)
        user_id = auth_client.get_user(cls.token)
        # WARNING: don't call any logging methods on the context object,
        # it'll result in a NoneType error
        cls.ctx = MethodContext(None)
        cls.ctx.update({'token': cls.token,
                        'user_id': user_id,
                        'provenance': [
                            {'service': 'kb_flye',
                             'method': 'please_never_use_it_in_production',
                             'method_params': []
                             }],
                        'authenticated': 1})
        cls.wsURL = cls.cfg['workspace-url']
        cls.wsClient = Workspace(cls.wsURL)
        cls.shockURL = cls.cfg['shock-url']
        cls.hs = HandleService(url=cls.cfg['handle-service-url'],
                               token=cls.token)
        cls.wsClient = Workspace(cls.wsURL)
        cls.serviceImpl = kb_flye(cls.cfg)
        cls.scratch = cls.cfg['scratch']
        cls.callback_url = os.environ['SDK_CALLBACK_URL']
        suffix = int(time.time() * 1000)
        cls.wsName = "test_flye_" + str(suffix)
        cls.wsinfo = cls.wsClient.create_workspace({'workspace': cls.wsName})
        print('created workspace ' + cls.getWsName())

        cls.readUtilsImpl = ReadsUtils(cls.callback_url, token=cls.token)
        cls.dfuClient = DataFileUtil(url=cls.callback_url, token=cls.token)
        cls.assemblyUtil = AssemblyUtil(cls.callback_url)
        cls.staged = {}
        cls.nodes_to_delete = []
        cls.handles_to_delete = []
        cls.setupTestData()
        print('\n\n=============== Starting Flye tests ==================')

    @classmethod
    def tearDownClass(cls):
        print('\n\n=============== Cleaning up ==================')
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')
        if hasattr(cls, 'nodes_to_delete'):
            for node in cls.nodes_to_delete:
                cls.delete_shock_node(node)
        if hasattr(cls, 'handles_to_delete'):
            if cls.handles_to_delete:
                cls.hs.delete_handles(cls.hs.hids_to_handles(cls.handles_to_delete))
                print('Deleted handles ' + str(cls.handles_to_delete))

    @classmethod
    def getWsName(cls):
        return cls.wsinfo[1]

    def getImpl(self):
        return self.serviceImpl

    @classmethod
    def delete_shock_node(cls, node_id):
        header = {'Authorization': 'Oauth {0}'.format(cls.token)}
        requests.delete(cls.shockURL + '/node/' + node_id, headers=header,
                        allow_redirects=True)
        print('Deleted shock node ' + node_id)

    @classmethod
    def upload_file_to_shock_and_get_handle(cls, test_file):
        '''                                                                                    
        Uploads the file in test_file to shock and returns the node and a                      
        handle to the node.                                                                    
        '''
        # file can't be in /kb/module/test or dfu won't find it                                
        temp_file = os.path.join("/kb/module/work/tmp", os.path.basename(test_file))
        shutil.copy(os.path.join("/kb/module/test", test_file), temp_file)

        print('loading file to shock: ' + test_file)
        fts = cls.dfuClient.file_to_shock({'file_path': temp_file,
                                           'make_handle':True})

        cls.nodes_to_delete.append(fts['shock_id'])
        cls.handles_to_delete.append(fts['handle']['hid'])

        return fts['shock_id'], fts['handle']['hid'], fts['size']

    @classmethod
    def upload_reads(cls, wsobjname, object_body, fwd_reads,
                     rev_reads=None, single_end=False, sequencing_tech='Illumina',
                     single_genome='1'):

        ob = dict(object_body)  # copy                                                         
        ob['sequencing_tech'] = sequencing_tech
        ob['wsname'] = cls.getWsName()
        ob['name'] = wsobjname
        if single_end or rev_reads:
            ob['interleaved'] = 0
        else:
            ob['interleaved'] = 1
        print('\n===============staging data for object ' + wsobjname +
              '================')
        print('uploading forward reads file ' + fwd_reads['file'])
        fwd_id, fwd_handle_id, fwd_size = \
            cls.upload_file_to_shock_and_get_handle(fwd_reads['file'])

        ob['fwd_id'] = fwd_id
        rev_id = None
        rev_handle_id = None
        if rev_reads:
            print('uploading reverse reads file ' + rev_reads['file'])
            rev_id, rev_handle_id, rev_size = \
                cls.upload_file_to_shock_and_get_handle(rev_reads['file'])
            ob['rev_id'] = rev_id
        obj_ref = cls.readUtilsImpl.upload_reads(ob)
        objdata = cls.wsClient.get_object_info_new({
            'objects': [{'ref': obj_ref['obj_ref']}]
            })[0]
        cls.staged[wsobjname] = {'info': objdata,
                                 'ref': cls.make_ref(objdata),
                                 'fwd_node_id': fwd_id,
                                 'rev_node_id': rev_id,
                                 'fwd_handle_id': fwd_handle_id,
                                 'rev_handle_id': rev_handle_id
                                 }

    @classmethod
    def setupTestData(cls):
        print('Shock url ' + cls.shockURL)
        # print('WS url ' + cls.wsClient.url)                                                       
        # print('Handle service url ' + cls.hs.url)                                                 
        print('staging data')

        reads = {'file': 'data/E_coli_PacBio_40x.fastq.gz',
                 'name': 'reads.fastq',
                 'type': 'fastq'}
        cls.upload_reads('reads', {'single_genome': 1}, reads,
                         single_end=True, sequencing_tech="Pacbio")

        print('Data staged.')

    @classmethod
    def make_ref(self, object_info):
        return str(object_info[6]) + '/' + str(object_info[0]) + \
            '/' + str(object_info[4])


    # NOTE: According to Python unittest naming rules test method names should start from 'test'. # noqa
    def run_flye(self,
                 output_contigset_name,
                 pacbio_raw_reads = None):

        params = {'workspace_name': self.getWsName(),
                  'output_contigset_name': output_contigset_name,
                  'pacbio_raw_reads': pacbio_raw_reads
                  }

        ret = self.serviceImpl.run_kb_flye(self.ctx, params)[0]
        self.assertReportOK(ret, output_contigset_name)

    def assertReportOK(self, ret_obj, assembly_name):
        """                                                                                         
        assertReportAssembly: given a report object, check the object existence                     
        """
        report = self.wsClient.get_objects2({
                        'objects': [{'ref': ret_obj['report_ref']}]})['data'][0]
        self.assertEqual('KBaseReport.Report', report['info'][2].split('-')[0])
        self.assertEqual(1, len(report['data']['objects_created']))
        self.assertEqual('Assembled contigs',
                         report['data']['objects_created'][0]['description'])
        print("**************Report Message*************\n")
        print(report['data']['text_message'])


    # Uncomment to skip this test                                                                   
    def test_pacbio(self):
        self.run_flye( 'output_contigset_name',
                       pacbio_raw_reads=self.staged['reads']['ref'])
