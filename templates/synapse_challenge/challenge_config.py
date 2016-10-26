##-----------------------------------------------------------------------------
##
## challenge specific code and configuration
##
##-----------------------------------------------------------------------------


## A Synapse project will hold the assetts for your challenge. Put its
## synapse ID here, for example
## CHALLENGE_SYN_ID = "syn1234567"
CHALLENGE_SYN_ID = "syn2813581"

## Name of your challenge, defaults to the name of the challenge's project
CHALLENGE_NAME = "SMC-Het"

## Synapse user IDs of the challenge admins who will be notified by email
## about errors in the scoring script
# Chris Bare, Kyle Ellrot, and Minjeong Ko
ADMIN_USER_IDS = [1421212, 297328, 3323319]

## Each question in your challenge should have an evaluation queue through
## which participants can submit their predictions or models. The queues
## should specify the challenge project as their content source. Queues
## can be created like so:
##   evaluation = syn.store(Evaluation(
##     name="My Challenge Q1",
##     description="Predict all the things!",
##     contentSource="syn1234567"))
## ...and found like this:
##   evaluations = list(syn.getEvaluationByContentSource('syn3375314'))
## Configuring them here as a list will save a round-trip to the server
## every time the script starts.
evaluation_queues = [
{u'contentSource': u'syn2813581',
  u'createdOn': u'2015-06-16T18:40:16.436Z',
  u'description': u'Queue for workflow submission for the SMC-Het Challenge',
  u'etag': u'161b8f99-1dbc-4e6f-a8f4-76400b2413cf',
  u'id': u'4487063',
  u'name': u'SMC-Het-Challenge-evaluation-queue',
  u'ownerId': u'1421212',
  u'status': u'OPEN',
  u'submissionInstructionsMessage': u'See: https://www.synapse.org/#!Synapse:syn2813581/wiki/',
  u'submissionReceiptMessage': u'Thanks for submitting your workflow'}]
evaluation_queue_by_id = {q['id']:q for q in evaluation_queues}

## define the default set of columns that will make up the leaderboard
LEADERBOARD_COLUMNS = [
    dict(name='objectId',      display_name='ID',      columnType='STRING', maximumSize=20),
    dict(name='userId',        display_name='User',    columnType='STRING', maximumSize=20, renderer='userid'),
    dict(name='entityId',      display_name='Entity',  columnType='STRING', maximumSize=20, renderer='synapseid'),
    dict(name='versionNumber', display_name='Version', columnType='INTEGER'),
    dict(name='name',          display_name='Name',    columnType='STRING', maximumSize=240),
    dict(name='team',          display_name='Team',    columnType='STRING', maximumSize=240)]

## Here we're adding columns for the output of our scoring functions, score,
## rmse and auc to the basic leaderboard information. In general, different
## questions would typically have different scoring metrics.
leaderboard_columns = {}
for q in evaluation_queues:
    leaderboard_columns[q['id']] = LEADERBOARD_COLUMNS + [
        dict(name='score',         display_name='Score',   columnType='DOUBLE'),
        dict(name='rmse',          display_name='RMSE',    columnType='DOUBLE'),
        dict(name='auc',           display_name='AUC',     columnType='DOUBLE')]

## map each evaluation queues to the synapse ID of a table object
## where the table holds a leaderboard for that question
leaderboard_tables = {}


def validate_submission(evaluation, submission):
    """
    Find the right validation function and validate the submission.

    :returns: (True, message) if validated, (False, message) if
              validation fails or throws exception
    """
    return True, "Looks OK to me!"


def score_submission(evaluation, submission):
    """
    Find the right scoring function and score the submission

    :returns: (score, message) where score is a dict of stats and message
              is text for display to user
    """
    import random
    return (dict(score=random.random(), rmse=random.random(), auc=random.random()), "You did fine!")


