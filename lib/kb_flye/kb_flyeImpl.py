# -*- coding: utf-8 -*-
#BEGIN_HEADER
from __future__ import print_function
import os
import re
import uuid
import requests
import json
import subprocess
import numpy as np
import yaml
import time
import zipfile
from pprint import pformat
import sys
from html import escape
from shutil import copy, copytree, move

from installed_clients.WorkspaceClient import Workspace
from installed_clients.ReadsUtilsClient import ReadsUtils  # @IgnorePep8
from installed_clients.baseclient import ServerError
from installed_clients.AssemblyUtilClient import AssemblyUtil
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.kb_quastClient import kb_quast
#END_HEADER


class kb_flye:
    '''
    Module Name:
    kb_flye

    Module Description:
    A KBase module: kb_flye
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "1.0.0"
    GIT_URL = ""
    GIT_COMMIT_HASH = ""

    #BEGIN_CLASS_HEADER
    def log(self, target, message):
        if target is not None:
            target.append(message)
        print(message)
        sys.stdout.flush()

    # from kb_SPAdes/utils/spades_utils.py:
    def load_stats(self, console, input_file_name):
        self.log(console, 'Starting conversion of FASTA to KBaseGenomeAnnotations.Assembly')
        # self.log(console, 'Building Object.')
        if not os.path.isfile(input_file_name):
            raise Exception('The input file name {0} is not a file!'.format(input_file_name))
        with open(input_file_name, 'r') as input_file_handle:
            contig_id = None
            sequence_len = 0
            length_dict = dict()
            coverage_dict = dict()
            first_header_found = False
            # Pattern for replacing white space
            pattern = re.compile(r'\s+')
            for current_line in input_file_handle:
                if (current_line[0] == '>'):
                    # found a header line
                    # Wrap up previous fasta sequence
                    if not first_header_found:
                        first_header_found = True
                    else:
                        length_dict[contig_id] = sequence_len
                        sequence_len = 0
                    fasta_header = current_line.replace('>', '').strip()
                    # self.log(console, 'fasta header = '+fasta_header)
                    try:
                        fields = fasta_header.strip().split(' ')
                        contig_id = fields[0]
                        # don't trust length from header, we look at seqence:
                        # sequence_len = int(fields[1][7:]) if (fields[1].startswith('length=')) else 0
                        coverage = float(
                            fields[2][6:-1]) if (fields[2].startswith('depth=')) else 0.0
                        # length_dict[contig_id] = sequence_len
                        coverage_dict[contig_id] = coverage
                    except (IndexError, ValueError, KeyError):
                        contig_id = fasta_header.strip()
                        coverage_dict[contig_id] = 0
                else:
                    sequence_len += len(re.sub(pattern, '', current_line))
        # wrap up last fasta sequence
        if not first_header_found:
            raise Exception("There are no contigs in this file")
        else:
            length_dict[contig_id] = sequence_len
        return [length_dict, coverage_dict]

    # from kb_SPAdes/utils/spades_utils.py:
    def mkdir_p(self, path):
        """
        mkdir_p: make directory for given path
        """
        if not path:
            return
        try:
            os.makedirs(path)
        except OSError as exc:
            if exc.errno == errno.EEXIST and os.path.isdir(path):
                pass
            else:
                raise

    def read_template(self, template_name):
        '''
        read in a template file and escape all html content
        used to display template contents
        '''
        with open(os.path.join(self.appdir, 'templates', template_name)) as file:
            lines = file.read()

        # escape all the html, display the results
        escaped_lines = escape(lines, quote=True)
        return escaped_lines

    def read_html(self, html_file):
        '''
        read in a html file
        '''
        with open(html_file) as file:
            lines = file.read()
        return lines

    # from kb_SPAdes/utils/spades_utils.py:
    def zip_folder(self, folder_path, output_path):
        """
        zip_folder: Zip the contents of an entire folder (with that folder included
        in the archive). Empty subfolders could be included in the archive as well
        if the commented portion is used.
        """
        with zipfile.ZipFile(output_path, 'w',
                             zipfile.ZIP_DEFLATED,
                             allowZip64=True) as ziph:
            for root, folders, files in os.walk(folder_path):
                for f in files:
                    absolute_path = os.path.join(root, f)
                    relative_path = os.path.join(os.path.basename(root), f)
                    # print "Adding {} to archive.".format(absolute_path)
                    ziph.write(absolute_path, relative_path)

        print("{} created successfully.".format(output_path))

    # from kb_SPAdes/utils/spades_utils.py:
    def generate_output_file_list(self, console, out_dir):
        """
        _generate_output_file_list: zip result files and generate file_links for report
        """
        self.log(console, 'start packing result files')

        output_files = list()

        output_directory = os.path.join(self.scratch, str(uuid.uuid4()))
        self.mkdir_p(output_directory)
        flye_output = os.path.join(output_directory, 'flye_output.zip')
        self.zip_folder(out_dir, flye_output)

        output_files.append({'path': flye_output,
                             'name': os.path.basename(flye_output),
                             'label': os.path.basename(flye_output),
                             'description': 'Output file(s) generated by Flye'})

        return output_files

    # adapted from kb_SPAdes/utils/spades_utils.py;
    # add templated report
    def generate_report(self, console, warnings, fa_file_name, params, out_dir, wsname):
        """
        Generating and saving report
        """
        self.log(console, 'Generating and saving report')

        fa_file_with_path = os.path.join(out_dir, fa_file_name)
        [length_stats, coverage_stats] = self.load_stats(console, fa_file_with_path)
        lengths = [length_stats[contig_id] for contig_id in length_stats]

        assembly_ref = wsname + '/' + params['output_contigset_name']

        report_text = ''
        report_text += 'Flye results saved to: ' + wsname + '/' + out_dir + '\n'
        report_text += 'Assembly saved to: ' + assembly_ref + '\n'
        report_text += 'Assembled into ' + str(len(lengths)) + ' contigs.\n'
        report_text += 'Avg Length: ' + str(sum(lengths) / float(len(lengths))) + ' bp.\n'

        # compute a simple contig length distribution
        bins = 10
        counts, edges = np.histogram(lengths, bins)
        report_text += 'Contig Length Distribution (# of contigs -- min to max ' + 'basepairs):\n'
        for c in range(bins):
            report_text += ('   ' + str(counts[c]) + '\t--\t' + str(edges[c]) + ' to ' +
                            str(edges[c + 1]) + ' bp\n')
        self.log(console, 'Running QUAST')
        kbq = kb_quast(self.callbackURL)
        quastret = kbq.run_QUAST(
            {'files': [{'path': fa_file_with_path, 'label': params['output_contigset_name']}]})
        # self.log(console,'quastret = '+pformat(quastret))

        # delete assembly file to keep it out of zip
        os.remove(fa_file_with_path)

        # make data table for report
        contig_data = []
        for contig_id in length_stats:
            contig_data.append({'contig_id': contig_id,
                                'coverage': coverage_stats[contig_id],
                                'length': length_stats[contig_id]})

        # self.log(console, 'contig_data = '+pformat(contig_data))

        # move quast output into main out_dir
        move(os.path.join(quastret['quast_path'], 'report.html'),
             os.path.join(out_dir, 'quast_report.html'))

        output_files = self.generate_output_file_list(console, out_dir)

        # render template
        template_file = 'flye_tabs.tt'
        tmpl_data = {
            'page_title': 'Flye Report',
            'data_array': contig_data,
            'cols': [
                {'data': 'contig_id',  'title': 'Contig ID'},
                {'data': 'coverage',   'title': 'Coverage (x)'},
                {'data': 'length',   'title': 'Length (bp)'}
            ]
        }
        # tmpl_data['quast_output'] = '<iframe>'+self.read_html(os.path.join(quastret['quast_path'],'report.html'))+'</iframe>'
        # tmpl_data['quast_output'] = '<iframe frameborder="0" width="100%" height="100%" src="'+os.path.join(quastret['quast_path'],'report.html')+'"></iframe>'
        tmpl_data['quast_output'] = '<iframe style="display:block; width:100%; height:100vh; border:none;" src="quast_report.html"></iframe>'
        tmpl_data['tmpl_vars'] = json.dumps(tmpl_data, sort_keys=True, indent=2)
        tmpl_data['template_content'] = self.read_template(template_file)
        tmpl_data['flye_log'] = '<p><pre>'+'<br>'.join(filter(lambda line: not (
            line.startswith('tput') or line.lstrip().startswith('0 / ')), console))+'</pre></p>'

        # save report
        self.log(console, 'Saving report')
        report_file = 'flye_report.html'

        # copy the templates into 'scratch', where they can be accessed by KBaseReport
        try:
            copytree(
                os.path.join(self.appdir, 'templates'),
                os.path.join(self.scratch, 'templates')
            )
        except Exception as e:
            self.log(console, 'Exception copying tree. '+str(e))

        reportClient = KBaseReport(self.callbackURL)
        template_output = reportClient.render_template({
            'template_file': os.path.join(self.scratch, 'templates', template_file),
            'template_data_json': json.dumps(tmpl_data),
            'output_file': os.path.join(out_dir, report_file)
        })

        report_output = reportClient.create_extended_report(
            {'message': report_text,
             'objects_created': [{'ref': assembly_ref, 'description': 'Assembled contigs'}],
             'direct_html_link_index': 0,
             'file_links': output_files,
             'html_links': [{'path': out_dir,
                             'name': report_file,
                             'label': 'Flye report',
                             'description': 'description of template report'
                             }
                            ],
             'warnings': warnings,
             'report_object_name': 'kb_flye_report_' + str(uuid.uuid4()),
             'workspace_name': params['workspace_name']})

        return report_output['name'], report_output['ref']

    # get long reads
    def download_long(self, console, warnings, token, wsname, long_reads_libraries, min_long_read_length):
        try:
            wsClient = Workspace(self.workspaceURL, token=token)
        except Exception as e:
            raise ValueError("unable to instantiate wsClient. "+str(e))

        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I,
         WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)  # object_info tuple

        total_read_length = 0
        long_reads_paths = []
        i = 0
        for lib in long_reads_libraries:
            try:
                i += 1
                obj_id = {'ref': lib if '/' in lib else (wsname + '/' + lib)}
                lib_obj_info = wsClient.get_object_info_new({'objects': [obj_id]})[0]
                lib_obj_type = lib_obj_info[TYPE_I]
                # remove trailing version
                lib_obj_type = re.sub('-[0-9]+\.[0-9]+$', "", lib_obj_type)
                lib_ref = str(lib_obj_info[WSID_I])+'/' + \
                    str(lib_obj_info[OBJID_I])+'/'+str(lib_obj_info[VERSION_I])

                ruClient = ReadsUtils(url=self.callbackURL, token=token)
                self.log(console, "Getting long reads (from reads library object).\n")
                result = ruClient.download_reads({'read_libraries': [lib_ref],
                                                  'interleaved': 'false'})
                long_reads_path = result['files'][lib_ref]['files']['fwd']
                [n_reads, n_reads_short, total_lib_read_length] = self.filter_short_fastq(
                    console, long_reads_path, min_long_read_length)
                total_read_length += total_lib_read_length
                    
                if (n_reads_short > 0):
                    self.log(warnings, "Warning:  Of "+str(n_reads)+" long reads, "+str(n_reads_short)+" are shorter than " +
                             str(min_long_read_length)+"; consider using the filtlong app to filter out shorter reads.")
                # work around minimap2 bug with long read names:
                long_reads_path = self.rename_fastq(console, long_reads_path, i)
                long_reads_paths.append(long_reads_path)
                    
            except Exception as e:
                raise ValueError('Unable to download long reads\n' + str(e))
        if (total_read_length > 10000000000):
            raise ValueError('Too many long reads; total length is limited to 10 GB and you have '+str(total_read_length)+' B.  Use filtlong app to filter out lower quality reads.')
        return long_reads_paths

    # examine fastq files, count total read length
    def filter_short_fastq(self, console, fastq_path, min_length):
        n_reads = 0
        n_reads_short = 0
        total_read_length = 0
        with open(fastq_path, 'r') as input_file_handle:
            for current_line in input_file_handle:
                if (current_line[0] == '@'):
                    # self.log(console, 'fastq header = '+current_line)
                    n_reads += 1
                    seq = next(input_file_handle)
                    if len(seq) < min_length:
                        n_reads_short += 1
                    total_read_length += len(seq)
                    next(input_file_handle)
                    next(input_file_handle)
            self.log(console, str(n_reads)+' long reads found, ' +
                     str(n_reads_short)+' under '+str(min_length)+' bp')
        return [n_reads, n_reads_short, total_read_length]

    # examine fastq files, count total read length
    def rename_fastq(self, console, fastq_path, index):
        directory, file_name = os.path.split(fastq_path)
        renamed_fastq_path = os.path.join(directory, 'renamed_'+file_name)
        cmd = 'seqtk rename '+fastq_path+' short'+str(index)+'- > '+renamed_fastq_path
        self.log(console, "command: "+cmd)
        cmdProcess = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                      stderr=subprocess.STDOUT, shell=True)
        for line in cmdProcess.stdout:
            self.log(console, line.decode("utf-8").rstrip())
        cmdProcess.wait()
        if cmdProcess.returncode != 0:
            raise ValueError('Error running '+cmd+'; see logs under "Job Status" for details.')
        return renamed_fastq_path

    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.cfg = config
        self.cfg['SDK_CALLBACK_URL'] = os.environ['SDK_CALLBACK_URL']
        self.cfg['KB_AUTH_TOKEN'] = os.environ['KB_AUTH_TOKEN']
        self.callbackURL = self.cfg['SDK_CALLBACK_URL']
        self.workspaceURL = config['workspace-url']
        self.shockURL = config['shock-url']
        self.scratch = os.path.abspath(config['scratch'])
        if not os.path.exists(self.scratch):
            os.makedirs(self.scratch)
        self.appdir = os.path.abspath(config['appdir'])
        #END_CONSTRUCTOR
        pass


    def run_kb_flye(self, ctx, params):
        """
        This example function accepts any number of parameters and returns results in a KBaseReport
        :param params: instance of mapping from String to unspecified object
        :returns: instance of type "ReportResults" -> structure: parameter
           "report_name" of String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN run_kb_flye
        console = []
        warnings = []
        self.log(console, 'Running run_flye with params:\n{}'.format(
            json.dumps(params, indent=1)))
        token = self.cfg['KB_AUTH_TOKEN']

        # param checks
        required_params = ['workspace_name',
                           'long_reads_libraries',
                           'long_reads_type',
                           'output_contigset_name']

        for required_param in required_params:
            if required_param not in params or params[required_param] is None:
                raise ValueError("Must define required param: '"+required_param+"'")

        if params['long_reads_type'] not in ["pacbio-raw", "pacbio-corr", "pacbio-hifi", "nano-raw", "nano-corr", "nano-hq"]:
            raise ValueError("long reads type '"+
                             str(params['long_reads_type'])+
                             "' not supported by Flye")

        if len(params['long_reads_libraries']) == 0:
            raise ValueError("Must provide at least one long reads library")

        # load provenance
        provenance = [{}]
        if 'provenance' in ctx:
            provenance = ctx['provenance']
        if 'input_ws_objects' not in provenance[0]:
            provenance[0]['input_ws_objects'] = []
        provenance[0]['input_ws_objects'].extend(params['long_reads_libraries'])

        # build command line
        cmd = '/kb/module/Flye-2.9.4/bin/flye'

        # download long libraries
        if 'long_reads_libraries' in params and params['long_reads_libraries'] is not None:
            longLibs = self.download_long(
                console, warnings, token, params['workspace_name'], params['long_reads_libraries'], 1000)
            cmd += ' --'+str(params['long_reads_type'])
            for longLib in longLibs:
                cmd += ' '+longLib

        if ('meta' in params and (params['meta'] == 1)):
            cmd += ' --meta'

        if 'min_overlap' in params and params['min_overlap'] is not None:
            cmd += ' --min-overlap '+str(params['min_overlap'])

        # output directory
        outputDir = os.path.join(self.scratch, "flye_"+str(uuid.uuid4()))
        self.mkdir_p(outputDir)
        cmd += ' --out-dir '+outputDir

        # run it
        self.log(console, "command: "+cmd)
        cmdProcess = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                      stderr=subprocess.STDOUT, shell=True)
        for line in cmdProcess.stdout:
            self.log(console, line.decode("utf-8").rstrip())
        cmdProcess.wait()
        if cmdProcess.returncode != 0:
            raise ValueError('Error running '+cmd+'; see logs under "Job Status" for details.')

        # save assembly
        try:
            contigsPath = os.path.join(outputDir, 'assembly.fasta')
            auClient = AssemblyUtil(url=self.callbackURL, token=token, service_ver='release')
            auClient.save_assembly_from_fasta(
                {'file': {'path': contigsPath},
                 'workspace_name': params['workspace_name'],
                 'assembly_name': params['output_contigset_name']})
        except Exception as e:
            raise ValueError('Error saving assembly\n' + str(e))

        # make report
        report_name, report_ref = self.generate_report(
            console, warnings, contigsPath, params, outputDir, params['workspace_name'])
        output = {'report_name': report_name,
                  'report_ref': report_ref}
        #END run_kb_flye

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method run_kb_flye return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
