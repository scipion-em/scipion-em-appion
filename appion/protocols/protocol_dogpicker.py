# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Laura del Cano (ldelcano@cnb.csic.es) [1]
# *              J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [2]
# *              Yunior C. Fonseca Reyna (cfonseca@cnb.csic.es) [1]
# *
# * [1] BCU, Centro Nacional de Biotecnologia, CSIC
# * [2] SciLifeLab, Stockholm University
# *
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import os

import pyworkflow.protocol.params as params
import pyworkflow.utils as pwutils
from pyworkflow.utils.properties import Message
from pwem.protocols import ProtParticlePickingAuto
from pwem.emlib.image import ImageHandler

from ..convert import readSetOfCoordinates


class DogPickerProtPicking(ProtParticlePickingAuto):
    """
    Automatically detects particle coordinates in cryo-EM micrographs using
    the Appion DogPicker approach. The protocol is designed to identify
    candidate particles by enhancing image features through Difference of
    Gaussian filtering, allowing rapid particle localization in large
    micrograph datasets.

    AI Generated:

    DogPicker Particle Picking (DogPickerProtPicking) - User Manual
        Overview

        The DogPicker Particle Picking protocol provides an automated method
        for identifying particles in cryo-EM micrographs using the Appion
        DogPicker framework. Its primary purpose is to accelerate particle
        selection during single-particle analysis workflows by detecting
        candidate molecular projections directly from raw or preprocessed
        micrographs.

        In practical cryo-EM processing pipelines, particle picking is one
        of the earliest and most influential stages because the quality of
        the selected coordinates strongly affects all downstream analyses,
        including 2D classification, 3D reconstruction, and refinement.
        This protocol is particularly useful for users who need a rapid and
        relatively simple automated picking strategy that can process large
        collections of micrographs consistently.

        Biological Motivation and Use Cases

        Biological cryo-EM datasets frequently contain thousands of
        particles distributed across many micrographs. Manual selection is
        often impractical, time-consuming, and subject to operator bias.
        Automated picking methods such as DogPicker allow users to generate
        reproducible particle coordinates suitable for initial processing
        and exploratory analysis.

        The protocol is commonly used during early stages of structure
        determination when users need a first estimate of particle
        populations. It can also serve as a preparatory step before more
        advanced deep learning or template-based picking strategies.
        Because the method relies primarily on local intensity features,
        it performs best when particles display relatively strong contrast
        against the surrounding background.

        Inputs and General Workflow

        The protocol requires a set of cryo-EM micrographs as input.
        These micrographs should ideally be motion-corrected and reasonably
        free of severe artifacts before particle detection begins. The
        quality of the input data strongly influences the reliability of
        the resulting particle coordinates.

        During processing, each micrograph is analyzed independently to
        identify local image features compatible with the expected particle
        appearance and size. The detected coordinates are then stored as a
        particle coordinate set that can be used directly in downstream
        extraction and classification workflows.

        Particle Diameter and Biological Interpretation

        One of the most important parameters is the expected particle
        diameter. This value determines the approximate spatial scale at
        which the protocol searches for particle-like features. Selecting
        an appropriate diameter is essential for biologically meaningful
        results.

        If the diameter is underestimated, large particles may be detected
        only partially or fragmented into multiple coordinates. If it is
        overestimated, neighboring particles or background contaminants may
        be merged incorrectly. In practice, the diameter should reflect the
        approximate maximum width of the molecular projection visible in
        the micrographs.

        For elongated or flexible complexes, users should generally choose
        a diameter corresponding to the dominant compact region rather than
        the longest dimension of the particle.

        Threshold Selection and Detection Sensitivity

        The detection threshold controls how strongly image features must
        stand out from the local background in order to be considered valid
        particles. Lower thresholds increase sensitivity and usually detect
        more candidate particles, but they also raise the number of false
        positives arising from contamination, carbon edges, or noise.

        Higher thresholds produce cleaner coordinate sets with fewer
        spurious detections, although genuine low-contrast particles may
        be missed. In biological practice, users often begin with moderate
        threshold values and then inspect the coordinates visually before
        deciding whether additional refinement is needed.

        The optimal threshold depends strongly on ice thickness, particle
        contrast, detector quality, and specimen heterogeneity.

        Image Contrast and Inversion

        The protocol optionally supports image inversion before particle
        detection. This option is biologically important because some
        datasets contain particles that appear darker than the background,
        while others display the opposite contrast relationship depending
        on imaging conditions and preprocessing conventions.

        Correctly matching the expected particle contrast improves picking
        reliability considerably. Users should visually inspect their
        micrographs before processing and determine whether particles
        appear brighter or darker relative to the surrounding vitreous ice.

        Additional Detection Parameters

        Advanced users may provide additional detection settings to adapt
        the protocol to unusually challenging datasets. These options can
        influence the number of detection scales explored, the allowable
        particle area, or the maximum number of candidate peaks identified
        in each micrograph.

        Such advanced adjustments are most useful in datasets with strong
        contamination, broad particle size variability, crowded fields, or
        highly heterogeneous particle distributions. For routine biological
        processing, however, the default parameters are often sufficient to
        obtain useful initial coordinate sets.

        Outputs and Their Interpretation

        After execution, the protocol generates a set of particle
        coordinates associated with the input micrographs. These
        coordinates define the estimated particle centers and are intended
        for downstream particle extraction.

        The resulting coordinate set should always be inspected visually.
        Automated particle picking methods can introduce both false
        positives and false negatives, especially in low-contrast or
        contaminated datasets. Biological interpretation should therefore
        rely on subsequent validation through 2D classification and
        refinement rather than on the picking results alone.

        Practical Recommendations

        In routine cryo-EM workflows, it is generally advisable to begin
        with conservative threshold values and a realistic particle
        diameter estimate. Visual inspection of the resulting coordinates
        on several representative micrographs provides the fastest way to
        determine whether the parameters are biologically appropriate.

        If too many contaminants are selected, increasing the threshold
        usually improves specificity. If genuine particles are missed,
        lowering the threshold moderately may improve recovery. Datasets
        containing strong ice gradients or carbon support regions often
        benefit from preprocessing or manual exclusion of problematic
        micrographs before automated picking.

        DogPicker is particularly effective as a fast initial picking tool
        for exploratory workflows and for generating particle sets that
        can later be cleaned through classification methods.

        Final Perspective

        Automated particle picking is a critical bridge between raw image
        acquisition and high-resolution structural analysis. Reliable
        coordinate detection allows cryo-EM users to process large
        datasets efficiently while reducing manual intervention. Careful
        selection of particle diameter, threshold sensitivity, and image
        contrast settings is essential for obtaining biologically useful
        particle populations suitable for downstream reconstruction and
        interpretation.
    """
    _label = 'dogpicker'

    def __init__(self, **args):
        ProtParticlePickingAuto.__init__(self, **args)

    # --------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):

        ProtParticlePickingAuto._defineParams(self, form)
        form.addParam('diameter', params.IntParam, default=100,
                      label='Diameter of particle in Å')
        form.addParam('invert', params.BooleanParam, default=False,
                      label='Invert',
                      help="Invert image before picking, DoG normally picks "
                           "white particles.")
        form.addParam('threshold', params.FloatParam, default=0.5,
                      label='Threshold',
                      help="Threshold in standard deviations above the mean, "
                           "e.g. --thresh=0.7")
        form.addParam('extraParams', params.StringParam,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Additional parameters',
                      help='Additional parameters for dogpicker: \n  '
                           '--num-slices=, --size-range=, --max-thresh=, --max-area='
                           '--max-peaks=')

    # --------------------------- STEPS functions ------------------------------
    def _pickMicrograph(self, mic, args):
        # Prepare mic folder and convert if needed
        micName = mic.getFileName()
        micDir = self._getTmpPath(pwutils.removeBaseExt(micName))
        pwutils.makePath(micDir)

        ih = ImageHandler()
        # If needed convert micrograph to mrc format, otherwise link it
        if pwutils.getExt(micName) != ".mrc":
            fnMicBase = pwutils.replaceBaseExt(micName, 'mrc')
            inputMic = os.path.join(micDir, fnMicBase)
            ih.convert(mic.getLocation(), inputMic)
        else:
            inputMic = os.path.join(micDir, os.path.basename(micName))
            pwutils.createLink(micName, inputMic)

        # Prepare environment
        from appion import Plugin
        Plugin.getEnviron()

        # Program to execute and it arguments
        program = "python2"
        outputFile = self._getExtraPath(pwutils.replaceBaseExt(inputMic, "txt"))

        args += " --image=%s --outfile=%s" % (inputMic, outputFile)

        dogpicker = Plugin.getHome("ApDogPicker.py")
        args = dogpicker + " " + args

        self.runJob(program, args)

    def createOutputStep(self):
        pass

    # --------------------------- INFO functions -------------------------------
    def _summary(self):
        summary = []
        summary.append("Number of input micrographs: %d"
                       % self.getInputMicrographs().getSize())
        if self.getOutputsSize() > 0:
            summary.append("Number of particles picked: %d"
                           % self.getCoords().getSize())
            summary.append("Particle size: %d" % self.getCoords().getBoxSize())
            summary.append("Threshold: %0.2f" % self.threshold)
            if self.extraParams.hasValue():
                summary.append("And other parameters: %s" % self.extraParams)
        else:
            summary.append(Message.TEXT_NO_OUTPUT_CO)
        return summary

    def _citations(self):
        return ['Voss2009']

    def _methods(self):
        methodsMsgs = []
        if self.getInputMicrographs() is None:
            return ['Input micrographs not available yet.']
        methodsMsgs.append("Input micrographs %s of size %d."
                           % (self.getObjectTag(self.getInputMicrographs()),
                              self.getInputMicrographs().getSize()))

        if self.getOutputsSize() > 0:
            output = self.getCoords()
            methodsMsgs.append('%s: User picked %d particles with a particle '
                               'size of %d and threshold %0.2f.'
                               % (self.getObjectTag(output), output.getSize(),
                                  output.getBoxSize(), self.threshold))
        else:
            methodsMsgs.append(Message.TEXT_NO_OUTPUT_CO)

        return methodsMsgs

    # --------------------------- UTILS functions ------------------------------
    def _getPickArgs(self):
        args = "--diam=%0.3f " % self.diameter.get()
        args += "--apix=%0.3f " % self.getInputMicrographs().getSamplingRate()
        args += "--thresh=%f" % self.threshold

        if self.invert:
            args += " --invert"

        args += " " + self.extraParams.get('')

        return [args]

    def readCoordsFromMics(self, workingDir, micList, coordSet):
        coordSet.setBoxSize(round(self.diameter.get() / self.getInputMicrographs().getSamplingRate()))
        readSetOfCoordinates(workingDir, micList, coordSet)

    def getCoordsDir(self):
        return self._getExtraPath()
