# -*- coding: utf-8 -*-
# Generated by Django 1.11.29 on 2024-04-17 11:52
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('pathogen_identification', '0043_rawreferencecompound'),
    ]

    operations = [
        migrations.AlterField(
            model_name='rawreferencecompound',
            name='family',
            field=models.ManyToManyField(blank=True, related_name='family', to='pathogen_identification.RawReference'),
        ),
        migrations.AlterField(
            model_name='rawreferencecompound',
            name='runs',
            field=models.ManyToManyField(blank=True, related_name='runs', to='pathogen_identification.RunMain'),
        ),
    ]